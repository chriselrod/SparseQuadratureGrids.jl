#Why the p in grid interfaces?
abstract type GridInterface{q <: QuadratureRule} end
abstract type CollapsedGrid{q <: QuadratureRule} end
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
struct FlattenedGrid{q} <: CollapsedGrid{q}
  nodes::Array{Float64,2}
  weights::Vector{Float64}
  baseline_weights::Vector{Float64}
  density::Vector{Float64}
end
struct FlatGrid{p,q} <: CollapsedGrid{q}
  nodes::Vector{SVector{p,Float64}}
  weights::Vector{Float64}
  baseline_weights::Vector{Float64}
  density::Vector{Float64}
end

struct SmolyakRaw{q, p, F} <: aPrioriBuild{q} end
struct Smolyak{q, p, F, T} <: aPrioriBuild{q} end

struct GridContainer{q} <: GridInterface{q}
  grids::Dict{Int, FlattenedGrid{q}}
  U::Array{Float64,2}
  l::Int
  seq::Vector{Int}
  mats::Dict{Int, Array{Float64,2}}
end
abstract type GridVessel{q, B <: GridBuild, K}  <: GridInterface{q} end
struct DynamicGridVessel{q, B, K} <: GridVessel{q, B, K}
  grids::Dict{K, FlatGrid{p,q} where p}
  U::Array{Float64, 2}
  mats::Dict{Int, Array{Float64,2}}
end
struct StaticGridVessel{q, B, K, p} <: GridVessel{q, B, K}
  grids::Dict{K, FlatGrid{p,q}}
  U::Array{Float64, 2}
  mats::Dict{Int, Array{Float64,2}}
end
function DynamicGridVessel(::Type{q},::Type{B},::Type{K},p::Int) where {q <: NestedQuadratureRule, B <: GridBuild, K}
  GridVessel{q,B,K}(Dict{K,FlatGrid{p,q} where p}(),Array{Float64,2}(p,p), Dict{Int,Array{Float64,2}}())
end
function StaticGridVessel(::Type{q},::Type{B},::Type{K},::Type{Val{p}}) where {q <: NestedQuadratureRule, B <: GridBuild, K, p}
  GridVessel{q,B,K,p}(Dict{K,FlatGrid{p,q}}(),Array{Float64,2}(p,p), Dict{Int,Array{Float64,2}}())
end



function aPrioriGrid(::Type{SmolyakRaw{q,p,F,T}}, f::Function, seq::Vector{Int}) where {q,p,F,T}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  eval_grid(FlatGrid{p,q}(grid), f)
end
function FlatGrid(::Type{SmolyakBuild{q,p,F}}, f::F, seq::Vector{Int}) where {q,p,F <: Function}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  eval_grid(FlatGrid{q}(grid), f)
end
function FlatGrid(::Type{SmolyakBuild{q,p,F}} where {F <: Function}, seq::Vector{Int}) where {p,q<:QuadratureRule}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  FlatGrid{q}(grid)
end
function eval_grid(g::FlatGrid{q}, f::Function)
  @inbounds for i ∈ eachindex(g.density)
    g.density[i] = g.baseline_weights[i] + f( g.nodes[:,i] )
  end
  normalize!(g)
end
function eval_grid(g::FlatGrid{q}, f::Function, ::Type{T})
  cache = Vector{T}(length(g.weights))
  @inbounds for i ∈ eachindex(g.density)
    f_val, cache[i] = f( g.nodes[:,i] )
    g.density[i] = g.baseline_weights[i] + f_val
  end
  normalize!(g).density, cache
end
function normalize!(g::FlatGrid)
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
  nodes = Array{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  density = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(ab.Grid.grid))
    weights[i] = ab.Grid.grid[j]
    density[i] = ab.F_cache[j] * weights[i]
    nodes[i] = get_node(ab, j)
  end
  density ./= sum(density)
  FlatGrid{q}(nodes, weights, Vector{Float64}(n), density)
end

function build(GV::GridVessel{q, B, K}, ::Type{Bc}, f::F, i::K, l::Int = 6) where {q,p,F<:Function, B <: AdaptiveBuild{q}, Bc <: AdaptiveBuild{q, p, F},K}
  ab = Adapt(Bc, f, l, i[end])
  GV.grids[i] = FlatGrid(ab)
  expand_grid!(ab.Grid, ab.Neighbors)
  density(ab)
end
function build(GV::GridVessel{q, B, K}, ::Type{Bc}, f::F, i::K) where {q,p,F<:Function, B <: aPrioriBuild{q}, Bc <: aPrioriBuild{q,p,F},K}
  GV.grids[i] = eval_grid(FlatGrid(Bc, i[end]), f)
end
function build(GV::GridVessel{q, B, K}, ::Type{Bc}, f::F, i::K) where {q,p,F<:Function, B <: aPrioriBuild{q}, Bc <: aPrioriBuild{q,p,F,T},K}
  grid = FlatGrid(Bc, i[end])
  GV.grids[i] = grid
  eval_grid(grid, f, T)
end

convert(::Type{FlattenedGrid{q}}, ab::AdaptiveBuild{q}) = FlattenedGrid(ab)
convert(::Type{FlatGrid{q}}, ab::AdaptiveBuild{q}) = FlatGrid(ab)

#calc_grid! performs a costly dynamic dispatch. The grids themselves are parametrized in their dimensionality.
#Benchmarks suggested it only takes a small handful of operations for the advantages of StaticArrays to outweigh the cost of a dynamic dispatch.
#It may be worthwhile to replace the [n x p] Array{Float64,2} with a [n] Vector{SVector{p}} in the FlattenedGrid.
#This would either increase the frequency of dynamic dispatces, or require some restructuring of the parametrization.
#If p=10, StaticArray flat grids could save > 400 microseconds (minus cost of extra dynamic dispatches).
#These dynamic dispatches could be avoided if we allow the option of forcing the Hessian to be full rank.
#This suggests two API modes of full rank vs dimension reduction.
#Full-rank forcing provides the advantage of preventing silent failures of the optimization step to find the posterior unconstrained mode.
#No additional dynamic dispatches needed on first call; only on future "eval_grid" calls.
function calc_grid!(GV::DynamicGridVessel{q, B, K} where {p, q <: QuadratureRule}, i::K, f::F ) where {B <: GridBuild, F <: Function, K}
  build!(GV, B{i[1],F}, f, i)#::(FlatGrid{p,q} where p) No point annotating, because the return remains dynamic.
end
function calc_grid!(GV::DynamicGridVessel{q, B, K} where {p, q <: QuadratureRule}, i::K, f::F, ::Type{T}) where {B <: GridBuild, F <: Function, T}
  build!(GV, B{i[1],F,T}, f, i)::Tuple{Vector{Float64},Vector{P}}
end

function calc_grid!(GV::StaticGridVessel{q, B, K, p} where {p, q <: QuadratureRule}, i::K, f::F ) where {B <: GridBuild, F <: Function}
  build!(GV, B{p,F}, f, i)::FlatGrid{p,q}
end
function calc_grid!(GV::StaticGridVessel{q, B, K, p} where {p, q <: QuadratureRule}, i::K, f::F, ::Type{T}) where {B <: GridBuild, F <: Function, T}
  build!(GV, B{p,F,T}, f, i)#::Tuple{Vector{Float64},Vector{P}} #Should be inferred correctly.
end





function FlattenedGrid(::Type{Val{p}}, l::Int, ::Type{q}, seq::Vector{Int} = default(q)) where {q}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, l)
  FlattenedGrid(grid)
end
function FlatGrid(::Type{Val{p}}, l::Int, ::Type{q}, seq::Vector{Int} = default(q)) where {q}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, l)
  FlatGrid(grid)
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
function flat_nodes_and_weights(Grid::NestedGrid{p,q})
  n = length(Grid.grid)
  nodes = Vector{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(Grid.grid))
    weights[i] = Grid.grid[j]
    nodes[i] = get_node.(Grid, j)
  end
  nodes, weights
end
squared_sum(x::AbstractArray) = sum(x .^ 2)
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
function FlatGrid(Grid::NestedGrid{p,q}) where {p,q<:NestedQuadratureRule}
  nodes, weights = flat_nodes_and_weights(Grid)
  Flat(nodes, weights, q)
end
FlatGrid(n::Vector{SVector{p,Float64}} where p, w::Vector{Float64}, ::Type{GenzKeister}) = FlatGrid(n, w, - squared_sum.(n), GenzKeister)
FlatGrid(n::Vector{SVector{p,Float64}} where p, w::Vector{Float64}, ::Type{KronrodPatterson}) = FlatGrid(n, w, zeros(w), KronrodPatterson)
function FlatGrid(n::Vector{SVector{p,Float64}} where p, w::Vector{Float64}, b::Vector{Float64}, ::Type{q}) where {q}
  FlatGrid{p,q}(n, w, b, similar(b))
end
function FlatGrid(ab::AdaptiveBuild{GenzKeister,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Vector{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  baseline = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(ab.Grid.grid))
    weights[i] = ab.Grid.grid[j]
    nodes[i] = get_node(ab, j)
    baseline[i] = ab.baseline_cache[j]
  end
  FlatGrid{p,GenzKeister}(nodes, weights, baseline, Vector{Float64}(n))
end
function FlatGrid(ab::AdaptiveBuild{KronrodPatterson,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Vector{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(ab.Grid.grid))
    weights[i] = ab.Grid.grid[j]
    nodes[i] = get_node(ab, j)
  end
  FlatGrid{p,KronrodPatterson}(nodes, weights, Vector{Float64}(n), Vector{Float64}(n))
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
