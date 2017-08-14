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
struct MatrixCache{p}
  d::Dict{Int,Array{Float64,2}}
end

struct GridContainer{p,q} <: GridInterface{q}
  grids::Dict{Int, FlattenedGrid{q}}
  U::Array{Float64,2}
  l::Int
  seq::Vector{Int}
  mats::MatrixCache{p}
end
abstract type GridVessel{q, B <: GridBuild, K, p}  <: GridInterface{q} end
struct DynamicGridVessel{q, B, K, p} <: GridVessel{q, B, K, p}
  grids::Dict{K, FlatGrid{dr,q} where dr}
  U::Array{Float64, 2}
  mats::MatrixCache{p}
end
struct StaticGridVessel{q, B, K, p} <: GridVessel{q, B, K, p}
  grids::Dict{K, FlatGrid{p,q}}
  U::Array{Float64, 2}
  mats::MatrixCache{p}
end
function DynamicGridVessel(::Type{q},::Type{B},::Type{K},::Type{Val{p}}) where {q <: NestedQuadratureRule, B <: GridBuild, K, p}
  DynamicGridVessel{q,B,K,p}(Dict{K,FlatGrid{dr,q} where dr}(),Array{Float64,2}(p,p), MatrixCache(Val{p}))
end
function StaticGridVessel(::Type{q},::Type{B},::Type{K},::Type{Val{p}}) where {q <: NestedQuadratureRule, B <: GridBuild, K, p}
  StaticGridVessel{q,B,K,p}(Dict{K,FlatGrid{p,q}}(),Array{Float64,2}(p,p), MatrixCache(Val{p}))
end

MatrixCache(::Type{Val{p}}) where p = MatrixCache{p}(Dict{Int,Array{Float64,2}}())
MatrixCache(p::Int) = MatrixCache{p}(Dict{Int,Array{Float64,2}}())

function FlatGrid(::Type{<:SmolyakBuild{q,p,F}}, f::F, seq::Vector{Int}) where {q,p,F <: Function}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  eval_grid(FlatGrid(grid), f)
end
function FlatGrid(::Type{<:SmolyakBuild{q,p,F}} where {F <: Function}, seq::Vector{Int}) where {p,q<:QuadratureRule}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, length(seq))
  FlatGrid(grid)
end
function eval_grid(g::FlatGrid, f::Function)
  @inbounds for i ∈ eachindex(g.density)
    g.density[i] = f( g.nodes[i] ) + g.baseline_weights[i]
  end
  normalize!(g)
end
function eval_grid(g::FlatGrid, f::Function, ::Type{T}) where T
  cache = Vector{T}( length( g.weights ) )
  @inbounds for i ∈ eachindex( g.density )
    f_val, cache[i] = f( g.nodes[i] )
    g.density[i] = f_val + g.baseline_weights[i]
  end
  normalize!(g).density, cache
end
function normalize!(g::FlatGrid)
  g.density .= exp.( g.density ) .* g.weights
  g.density ./= sum( g.density )
  g
end
function Adapt(::Type{Adaptive{q, p, F, T}}, f::F, l::Int, n::Int) where {q,p, F <: Function, T}
  ab = Adaptive{q, p, F, T}(Dict{SVector{p, Int},T}(), Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int}, Vector{Float64}}(), Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int},Float64}(), Dict{SVector{p, Int64}, Δprod{p}}(), NestedGrid(Val{p}, q), initialize_e(Val{p}), f, l)
  construct!(ab, n)
  ab
end
function Adapt(::Type{AdaptiveRaw{q, p, F}}, f::F, l::Int, n::Int) where {q,p, F <: Function}
  ab = AdaptiveRaw{q, p, F}(Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int}, Vector{Float64}}(), Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int},Float64}(), Dict{SVector{p, Int64}, Δprod{p}}(), NestedGrid(Val{p}, q), initialize_e(Val{p}), f, l)
  construct!(ab, n)
  ab
end
#Caution: this method is a dynamic dispatch.
Adaptive(p::Int, ::F, ::Type{T}, n::Int = 1_000, q::QuadratureRule = GenzKeister, l::Int = 6) where {F <: Function, T} = Adaptive(Adaptive{q, p, F, T}, l, n)

function density(ab::Adaptive{q, p, F, T}) where {q,p,F,T}
  n = length(ab.Grid.grid)
  density = Vector{Float64}(n)
  params = Vector{T}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    density[i] = ab.F_cache[j] * k
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
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    density[i] = ab.F_cache[j] * k
    nodes[i] = get_node(ab, j)
  end
  density ./= sum(density)
  FlatGrid{p,q}(nodes, weights, Vector{Float64}(n), density)
end

function build!(GV::GridVessel{q, B, K, fr} where fr, ::Type{Bc}, f::F, i::K, l::Int = 6) where {q, p, F<:Function, B <: AdaptiveBuild{q}, Bc <: AdaptiveBuild{q, p, F}, K}
  ab = Adapt(Bc, f, l, i[end])
  GV.grids[i] = FlatGrid(ab)
  for v ∈ values(ab.Neighbors)
    expand_grid!(ab.Grid, v)
  end
  density(ab)
end
function build!(GV::GridVessel{q, B, K, fr} where fr, ::Type{Bc}, f::F, i::K) where {q, p, F<:Function, B <: aPrioriBuild{q}, Bc <: aPrioriBuild{q,p,F}, K}
  GV.grids[i] = eval_grid(FlatGrid(Bc, ind(i)), f)
end
function build!(GV::GridVessel{q, B, K, fr} where fr, ::Type{Bc}, f::F, i::K) where {q, p, T, F<:Function, B <: aPrioriBuild{q}, Bc <: Smolyak{q,p,F,T}, K}
  grid = FlatGrid(Bc, ind(i))
  GV.grids[i] = grid
  eval_grid(grid, f, T)
end
ind(i::Vector{Int}) = i
ind(i::Tuple{Int,Vector{Int}}) = i[2]

convert(::Type{FlattenedGrid}, ab::AdaptiveBuild) = FlattenedGrid(ab)
convert(::Type{FlatGrid}, ab::AdaptiveBuild) = FlatGrid(ab)

#calc_grid! performs a costly dynamic dispatch. The grids themselves are parametrized in their dimensionality.
#Benchmarks suggested it only takes a small handful of operations for the advantages of StaticArrays to outweigh the cost of a dynamic dispatch.
#It may be worthwhile to replace the [n x p] Array{Float64,2} with a [n] Vector{SVector{p}} in the FlattenedGrid.
#This would either increase the frequency of dynamic dispatces, or require some restructuring of the parametrization.
#If p=10, StaticArray flat grids could save > 400 microseconds (minus cost of extra dynamic dispatches).
#These dynamic dispatches could be avoided if we allow the option of forcing the Hessian to be full rank.
#This suggests two API modes of full rank vs dimension reduction.
#Full-rank forcing provides the advantage of preventing silent failures of the optimization step to find the posterior unconstrained mode.
#No additional dynamic dispatches needed on first call; only on future "eval_grid" calls.
function calc_grid!(GV::DynamicGridVessel{q, B, K, fr} where {p, q <: QuadratureRule, fr}, i::K, f::F ) where {B <: GridBuild, F <: Function, K}
  build!(GV, B{i[1],F}, f, i)#::(FlatGrid{p,q} where p) No point annotating, because the return remains dynamic.
end
function calc_grid!(GV::DynamicGridVessel{q, B, K, fr} where {p, q <: QuadratureRule, fr}, i::K, f::F, ::Type{T}) where {B <: GridBuild, F <: Function, T, K}
  build!(GV, B{i[1],F,T}, f, i)::Tuple{Vector{Float64},Vector{P}}
end

function calc_grid!(GV::StaticGridVessel{q, B, K, p} where {q <: QuadratureRule}, i::K, f::F ) where {p, B <: GridBuild, F <: Function, K}
  build!(GV, B{p,F}, f, i)#::FlatGrid{p,q} #Should be inferred correctly.
end
function calc_grid!(GV::StaticGridVessel{q, B, K, p} where {q <: QuadratureRule}, i::K, f::F, ::Type{T}) where {p, B <: GridBuild, F <: Function, T, K}
  build!(GV, B{p,F,T}, f, i)#::Tuple{Vector{Float64},Vector{P}} #Should be inferred correctly.
end





function FlattenedGrid(::Type{Val{p}}, l::Int, ::Type{q}, seq::Vector{Int} = default(q)) where {p, q}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, l)
  FlattenedGrid(grid)
end
function FlatGrid(::Type{Val{p}}, l::Int, ::Type{q}, seq::Vector{Int} = default(q)) where {p, q}
  grid = NestedGrid(Val{p}, q, seq)
  smolyak!(grid, l)
  FlatGrid(grid)
end

function flatten_nodes_and_weights(Grid::NestedGrid{p,q}) where {p,q}
  n = length(Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(Grid.grid.w)
    weights[i] = k
    nodes[:,i] .= get_node.(Grid, j)
  end
  nodes, weights
end
function flat_nodes_and_weights(Grid::NestedGrid{p,q}) where {p,q}
  n = length(Grid.grid)
  nodes = Vector{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  for (i, (j,k)) ∈ enumerate(Grid.grid.w)
    weights[i] = k
    nodes[i] = get_node.(Grid, j)
  end
  nodes, weights
end

function FlattenedGrid(Grid::NestedGrid{p,q}) where {p,q<:NestedQuadratureRule}
  nodes, weights = flatten_nodes_and_weights(Grid)
  FlattenedGrid(nodes, weights, q)
end
FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, ::Type{GenzKeister}) = FlattenedGrid(n, w, vec(sum(abs2, n, [1])), GenzKeister)
FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, ::Type{KronrodPatterson}) = FlattenedGrid(n, w, zeros(w), KronrodPatterson)
function FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, b::Vector{Float64}, ::Type{q}) where {q}
  FlattenedGrid{q}(n, w, b, similar(b))
end
function FlattenedGrid(ab::AdaptiveBuild{GenzKeister,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  baseline = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    nodes[:,i] .= get_node(ab, j)
    baseline[i] = ab.baseline_cache[j]
  end
  FlattenedGrid{GenzKeister}(nodes, weights, baseline, Vector{Float64}(n))
end
function FlattenedGrid(ab::AdaptiveBuild{KronrodPatterson,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    nodes[:,i] .= get_node(ab, j)
  end
  FlattenedGrid{KronrodPatterson}(nodes, weights, Vector{Float64}(n), Vector{Float64}(n))
end
function FlatGrid(Grid::NestedGrid{p,q}) where {p,q<:NestedQuadratureRule}
  nodes, weights = flat_nodes_and_weights(Grid)
  FlatGrid(nodes, weights, q)
end
FlatGrid(n::Vector{SVector{p,Float64}} where p, w::Vector{Float64}, ::Type{GenzKeister}) = FlatGrid(n, w, sum.(abs2, n), GenzKeister)
FlatGrid(n::Vector{SVector{p,Float64}} where p, w::Vector{Float64}, ::Type{KronrodPatterson}) = FlatGrid(n, w, zeros(w), KronrodPatterson)
function FlatGrid(n::Vector{SVector{p,Float64}}, w::Vector{Float64}, b::Vector{Float64}, ::Type{q}) where {p,q}
  FlatGrid{p,q}(n, w, b, similar(b))
end
function FlatGrid(ab::AdaptiveBuild{GenzKeister,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Vector{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  baseline = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    nodes[i] = get_node(ab, j)
    baseline[i] = ab.baseline_cache[j]
  end
  FlatGrid{p,GenzKeister}(nodes, weights, baseline, Vector{Float64}(n))
end
function FlatGrid(ab::AdaptiveBuild{KronrodPatterson,p,F}) where {p,F}
  n = length(ab.Grid.grid)
  nodes = Vector{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    nodes[i] = get_node(ab, j)
  end
  FlatGrid{p,KronrodPatterson}(nodes, weights, Vector{Float64}(n), Vector{Float64}(n))
end


function getindex(g::GridContainer{p, q} where {p}, i::Int) where {q}
  get!(() -> FlattenedGrid(i, g.l, q, g.seq), g.grids, i)
end
mats(mc::MatrixCache{p}, r) where p = get!(() -> Array{Float64,2}(p,r), mc.d, r)
mats(g::GridVessel, r::Int) = mats(g.mats, r)

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
  @inbounds for (i, w) ∈ Grid.grid.w
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
  @views Grid.baseline_weights_positive .= sum(abs2, Grid.nodes_positive, [1])[1,:]
  @views Grid.baseline_weights_negative .= sum(abs2, Grid.nodes_negative, [1])[1,:]
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end

function set_baseline!(Grid::SplitWeights{p, KronrodPatterson}) where {p}
  fill!(Grid.baseline_weights_positive, 0.0)
  fill!(Grid.baseline_weights_negative, 0.0)
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end
