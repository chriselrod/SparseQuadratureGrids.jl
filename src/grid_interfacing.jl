#Why the p in grid interfaces?
abstract type GridInterface{q <: QuadratureRule} end
struct SplitWeights{p,q<: QuadratureRule}
  nodes_positive::Array{Float64,2}
  nodes_negative::Array{Float64,2}
  weights_positive::Array{Float64,1}
  weights_negative::Array{Float64,1}
  baseline_weights_positive::Array{Float64,1}
  baseline_weights_negative::Array{Float64,1}
  weight_sum::Array{Float64,1}
  total_count::Int64
end
struct FlattenedGrid{q <: QuadratureRule}
  nodes::Array{Float64,2}
  weights::Vector{Float64}
  baseline_weights::Vector{Float64}
  density::Vector{Float64}
end
#Maybe leave Smolyak abstract?
#GridFit only for aPriori types
#AdaptiveGrid is what you get for the posteriori builds.
struct GridFit{B <: aPrioriBuild}


end
struct SmolyakRaw{q, p, F} <: aPrioriBuild{q}
  Grid::NestedGrid{p, q}
  seq::Vector{Int}
  l::Int
  f::Function
end
struct Smolyak{q, p, F, T} <: aPrioriBuild{q}
  cache::Vector{T}
  Grid::NestedGrid{p, q}
  seq::Vector{Int}
  l::Int
  f::Function
end

struct GridContainer{q} <: GridInterface{q}
  grids::Dict{Int, FlattenedGrid{q}}
  U::Array{Float64,2}
  l::Int
  seq::Vector{Int}
  mats::Dict{Int, Array{Float64,2}}
end
struct GridVessel{q, B <: GridBuild} <: GridInterface{q}
  grids::Dict{K, FlattenedGrid{q}}
  U::Array{Float64, 2}
  mats::Dict{Int, Array{Float64,2}}
end
function GridVessel(::Type{q},::Type{B},::Type{K},p::Int) where {q <: NestedQuadratureRule, B <: GridBuild}
  GridVessel{q,B}(Dict{K,FlattenedGrid{q}}(),Array{Float64,2}(p,p) ,Dict{Int,Array{Float64,2}}())
end


function SmolyakRaw(::Type{Val{p}}, F::Function, q::QuadratureRule = GenzKeister, l::Int = 6, seq::Vector{Int} = default(q))
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, l)
  g = FlattenedGrid(grid)

  @inbounds for i ∈ eachindex(g.density)
    g.density[i] = g.baseline_weights[i] + F( g.nodes[:,i] )
  end

  g.density .= exp.(minimum(g.density) .- g.density) .* g.weights
  g.density ./= sum(g.density)

  SmolyakRaw{q, p, F}(g, seq, l, F)
end
function Smolyak(::Type{Val{p}}, F::Function, ::Type{T}, q::QuadratureRule = GenzKeister, l::Int = 6, seq::Vector{Int} = default(q))
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, l)
  Smolyak(FlattenedGrid(grid), F, T, l, seq)
end

function Smolyak(g::FlattenedGrid{q}, F::Function, ::Type{T}, l::Int = 6, seq::Vector{Int} = default(q))
  cache = Vector{T}(length(g.weights))

  @inbounds for i ∈ eachindex(g.density)
    f_val, cache[i] = F( g.nodes[:,i] )
    g.density[i] = g.baseline_weights[i] + f_val
  end
  normalize!(g)
  Smolyak{q, p, F}(cache, g, seq, l, F)
end

function aPrioriGrid(::Type{SmolyakRaw{q,p,F,T}}, f::Function, seq::Vector{Int}) where {q,p,F,T}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  eval_grid(FlattenedGrid{q}(grid), f)
end
function FlattenedGrid(::Type{SmolyakBuild{q,p,F}}, f::F, seq::Vector{Int}) where {F <: Function}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  eval_grid(FlattenedGrid{q}(grid), f)
end
function FlattenedGrid(::Type{SmolyakBuild{q,p,F}} where {q<:QuadratureRule, p, F <: Function}, seq::Vector{Int})
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  FlattenedGrid{q}(grid)
end
function eval_grid(g::FlattenedGrid{q}, f::Function)
  @inbounds for i ∈ eachindex(g.density)
    g.density[i] = g.baseline_weights[i] + f( g.nodes[:,i] )
  end
  g.density .= exp.(minimum(g.density) .- g.density) .* g.weights
  g.density ./= sum(g.density)
  normalize!(g)
end
function eval_grid(g::FlattenedGrid{q}, f::Function, ::Type{T})
  cache = Vector{T}(length(g.weights))
  @inbounds for i ∈ eachindex(g.density)
    f_val, cache[i] = f( g.nodes[:,i] )
    g.density[i] = g.baseline_weights[i] + f_val
  end
  normalize!(g).density, cache
end
function normalize!(g::FlattenedGrid)
  g.density .= exp.(minimum(g.density) .- g.density) .* g.weights
  g.density ./= sum(g.density)
  g
end
function Adapt(::Type{Adaptive{q, p, F, T}}, f::F, l::Int, n::Int) where {q,p, F <: Function, T}
  ab = Adaptive{q, p, F, T}(Dict{SVector{p, Int},T}(), Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int}, Vector{Float64}}(), Dict{SVector{p,Int}, Float64}(), Dict{SArray{p,Int},Float64}(), Dict{SVector{p, Int64}, Δprod{p}}(), NestedGrid{p,q}(Val{p}, q), initialize_e(p), f, l)
  construct!(ab, n)
  ab
end
function Adapt(::Type{AdaptiveRaw{q, p, F}}, f::F, l::Int, n::Int) where {q,p, F <: Function, T}
  ab = AdaptiveRaw{q, p, F}(Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int}, Vector{Float64}}(), Dict{SVector{p,Int}, Float64}(), Dict{SArray{p,Int},Float64}(), Dict{SVector{p, Int64}, Δprod{p}}(), NestedGrid(p, q), initialize_e(p), f, l)
  construct!(ab, n)
  ab
end
#Caution: this method is a dynamic dispatch.
Adaptive(p::Int, ::F, ::Type{T}, n::Int = 10_000, q::QuadratureRule = GenzKeister, l::Int = 6) where {F <: Function, T} = Adaptive(Adaptive{q, p, F, T}, l, n)

function density(ab::Adaptive{q, p, F, T}) where {q,p,F,T}
  n = length(ab.Grid.grid)
  density = Vector{Float64}(n)
  params = Vector{T}(n)
  for (i, j) ∈ enumerate(keys(ab.Grid.grid))
    density[i] = ab.F_cache[j] * ab.Grid.grid[j]
    params[i] = ab.cache[j]
  end
  density ./= sum(density)
  density, params
end
function density(ab::AdaptiveRaw{q, p, F}) where {q,p,F}
  n = length(ab.Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  density = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(ab.Grid.grid))
    weights[i] = ab.Grid.grid[j]
    density[i] = ab.F_cache[j] * weights[i]
    nodes[:,i] .= get_node(ab, j)
  end
  density ./= sum(density)
  FlattenedGrid{q}(nodes, weights, Vector{Float64}(n), density)
end

function build(GV::GridVessel{q, B}, ::Type{Bc}, f::F, i::Tuple{Int, Int, Int}, l::Int = 6) where {q,p,F<:Function, B <: AdaptiveBuild{q}, Bc <: AdaptiveBuild{q, p, F}}
  ab = Adapt(Bc, f, l, i[3])
  GV.grids[i] = FlattenedGrid(ab)
  expand_grid!(ab.Grid, ab.Neighbors)
  density(ab)
end
function build(GV::GridVessel{q, B}, ::Type{Bc}, f::F, i::Tuple{Int, Vector{Int}}) where {q,p,F<:Function, B <: aPrioriBuild{q}, Bc <: aPrioriBuild{q,p,F}}
  grid = FlattenedGrid(Bc, i[2])
  GV.grids[i] = grid
  eval_grid(grid, f)
end
function build(GV::GridVessel{q, B}, ::Type{Bc}, f::F, i::Tuple{Int, Vector{Int}}) where {q,p,F<:Function, B <: aPrioriBuild{q}, Bc <: aPrioriBuild{q,p,F,T}}
  grid = FlattenedGrid(Bc, i[2])
  GV.grids[i] = grid
  eval_grid(grid, f, T)
end

convert(::Type{FlattenedGrid{q}}, ab::AdaptiveBuild{q}) = FlattenedGrid(ab)

#calc_grid! performs a costly dynamic dispatch.
function calc_grid!(GV::GridVessel{p, q, B} where {p, q <: QuadratureRule}, i::Tuple{Int, Int}, f::F ) where {B <: GridBuild, F <: Function}
  build!(GV, B{i[1],F}, f, i)::FlattenedGrid{q}
end
function calc_grid!(GV::GridVessel{p, q, B} where {p, q <: QuadratureRule}, i::Tuple{Int, Int, Int}, f::F, ::Type{T}) where {B <: GridBuild, F <: Function, T}
  build!(GV, B{i[1],F,T}, f, i)::Tuple{Vector{Float64},Vector{P}}
end

###These are the two functions called by the JointPosteriors package.
###The first, return_grid!, is for the raw version, and it simply returns a flattened grid of the raw unconstrained values.
###The second returns a vector of parameter objects. More costly, but should be cheaper to compute marginals on thanks to cacheing the transformations.
function return_grid!(M::Model{q, P, B} where {q <: QuadratureRule, P <: parameters}, i::Tuple{Int, Vector{Int}}, f::F ) where {B <: GridBuild, F <: Function}
  haskey(GV.grids, i) ? eval_grid(GV.grids[i], f) : calc_grid!(GV, i, f)
end
function get!(M::Model{q, P, B} where {q <: QuadratureRule}, i::Tuple{Int, Int, Int}, f::F ) where {B <: GridBuild, F <: Function, P <: parameters}
  haskey(GV.grids, i) ? eval_grid(GV.grids[i], f, P) : calc_grid!(GV, i, f, P)
end



function FlattenedGrid(::Type{Val{p}}, l::Int, ::Type{q}, seq::Vector{Int} = default(q)) where {q}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, l)
  FlattenedGrid(grid)
end

function flatten_nodes_and_weights(Grid::NestedGrid{p,q})
  n = length(Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(Grid.grid))
    weights[i] = Grid.grid[j]
    nodes[:,i] .= get_node.(Grid, j)
  end
  nodes, weights
end

function FlattenedGrid(Grid::NestedGrid{p,q}) where {p,q<:NestedQuadratureRule}
  nodes, weights = flatten_nodes_and_weights(Grid)
  FlattenedGrid(nodes, weights, q)
end
FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, ::Type{GenzKeister}) = FlattenedGrid(n, w, - vec(sum(n .^ 2, [1])), GenzKeister)
FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, ::Type{KronrodPatterson}) = FlattenedGrid(n, w, zeros(w), KronrodPatterson)
function FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, b::Vector{Float64}, ::Type{q}) where {q}
  FlattenedGrid{q}(n, w, b, similar(b))
end
function FlattenedGrid(ab::AdaptiveBuild{GenzKeister,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  baseline = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(ab.Grid.grid))
    weights[i] = ab.Grid.grid[j]
    nodes[:,i] .= get_node(ab, j)
    baseline[i] = ab.baseline_cache[j]
  end
  FlattenedGrid{GenzKeister}(nodes, weights, baseline, Vector{Float64}(n))
end
function FlattenedGrid(ab::AdaptiveBuild{KronrodPatterson,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(ab.Grid.grid))
    weights[i] = ab.Grid.grid[j]
    nodes[:,i] .= get_node(ab, j)
  end
  FlattenedGrid{KronrodPatterson}(nodes, weights, Vector{Float64}(n), Vector{Float64}(n))
end


function getindex(g::GridContainer{p, q} where {p}, i::Int) where {q}
  get!(() -> FlattenedGrid(i, g.l, q, g.seq), g.grids, i)
end
function mats(g::GridInterface{p,q} where {q}, r::Int) where {p}
  get!(() -> Array{Float64,2}(p,r), g.mats, r)
end
GridContainer(p::Int,l::Int = 6,::Type{q} = GenzKeister,seq::Vector{Int} = default(q)) where {q} = GridContainer{p,q}(Dict{Tuple{Int,Int},FlattenedGrid{q}}(), zeros(p,p), l, seq, Dict{Int,Array{Float64,2}}())

function SplitWeights(p::Int64, total_count::Int64, negative_count::Int64, q::DataType)
  positive_count = total_count - negative_count
  SplitWeights{p,q}(Array{Float64,2}(p, positive_count), Array{Float64,2}(p, negative_count), Array{Float64,1}(positive_count), Array{Float64,1}(negative_count), Array{Float64,1}(positive_count), Array{Float64,1}(negative_count), Array{Float64,1}(2), total_count)
end

function SplitWeights(Grid::NestedGrid{p,q}) where {p, q <: NestedQuadratureRule}
#  nodes = maximum(values(Grid.Qs)).n
#  complete_grid = sum(values(Grid.Δ_prods))
  G = SplitWeights(p, length(Grid.grid), sum([i < 0 for i ∈ values(Grid.grid)]), q)
  nc = 0; pc = 0
  for i ∈ keys(Grid.grid)
    w = Grid.grid[i]
    w < 0 ? begin
      nc += 1
      G.nodes_negative[:,nc] .= get_node.(Grid, i)
      G.weights_negative[nc] = w
    end : begin
      pc += 1
      G.nodes_positive[:,pc] .= get_node.(Grid, i)
      G.weights_positive[pc] = w
    end
  end
  set_baseline!(G)
  G
end

function SplitWeights(p::Int, l::Int, ::Type{q}) where {q <: QuadratureRule}
  grid = NestedGrid(p, q)
  smolyak!(grid, l)
  SplitWeights(grid)
end
function SplitWeights(p::Int, l::Int, ::Type{q}, seq::Vector{Int}) where {q <: QuadratureRule}
  grid = NestedGrid(p, q, seq)
  smolyak!(grid, l)
  SplitWeights(grid)
end



function set_baseline!(Grid::SplitWeights{p, GenzKeister}) where {p}
  @views Grid.baseline_weights_positive .= - sum(Grid.nodes_positive .^ 2, [1])[1,:]
  @views Grid.baseline_weights_negative .= - sum(Grid.nodes_negative .^ 2, [1])[1,:]
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end

function set_baseline!(Grid::SplitWeights{p, KronrodPatterson}) where {p}
  fill!(Grid.baseline_weights_positive, 0.0)
  fill!(Grid.baseline_weights_negative, 0.0)
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end
