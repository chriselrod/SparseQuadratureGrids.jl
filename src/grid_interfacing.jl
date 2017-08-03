abstract type GridInterface{p,q <: QuadratureRule} end
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
struct SmolyakConstrained{q, p, F} <: aPrioriBuild{q}
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
struct GridContainer{p,q} <: GridInterface{p,q}
  grids::Dict{Int, FlattenedGrid{q}}
  U::Array{Float64,2}
  l::Int
  seq::Vector{Int}
  mats::Dict{Int,Array{Float64,2}}
end
struct GridVessel{p, q, B <: GridBuild} <: GridInterface{p,q}
  grids::Dict{Tuple{Int,Int}, FlattenedGrid{q}}
  U::Array{Float64, 2}
  seq::Vector{Int}
  mats::Dict{Int, Array{Float64,2}}
end

function AdaptiveConstrained(p::Int, q::QuadratureRule, F::Function, l::Int = 6)
  AdaptiveBuilder{p,q,F}(AdaptiveBuildInit(p, q, F, l)...)
end

function AdaptiveConstrained(p::Int, F::Function, n::Int = 10_000, q::QuadratureRule = GenzKeister, l::Int = 6)
  ab = AdaptiveConstrained(p, q, F, l)
  construct!(ab, n)
  ab
end


function SmolyakConstrained(p::Int, F::Function, q::QuadratureRule = GenzKeister, l::Int = 6, seq::Vector{Int} = default(q))
  grid = NestedGrid(p, q, seq)
  smolyak!(grid, l)
  g = FlattenedGrid(grid)

  @inbounds for i ∈ eachindex(g.density)
    g.density[i] = g.baseline_weights[i] + F( g.nodes[:,i] )
  end

  g.density .= exp.(minimum(g.density) .- g.density) .* g.weights
  g.density ./= sum(g.density)

  SmolyakConstrained{q, p, F}(g, seq, l, F)
end
function Smolyak(p::Int, F::Function, ::Type{T}, q::QuadratureRule = GenzKeister, l::Int = 6, seq::Vector{Int} = default(q))
  grid = NestedGrid(p, q, seq)
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
function eval_grid(g::FlattenedGrid{q}, F::Function)
  @inbounds for i ∈ eachindex(g.density)
    g.density[i] = g.baseline_weights[i] + F( g.nodes[:,i] )
  end
  g.density .= exp.(minimum(g.density) .- g.density) .* g.weights
  g.density ./= sum(g.density)
  normalize!(g)
end
function eval_grid_cache(g::FlattenedGrid{q}, F::Function, ::Type{T})
  cache = Vector{T}(length(g.weights))
  @inbounds for i ∈ eachindex(g.density)
    f_val, cache[i] = F( g.nodes[:,i] )
    g.density[i] = g.baseline_weights[i] + f_val
  end
  normalize!(g), cache
end
function normalize!(g::FlattenedGrid)
  g.density .= exp.(minimum(g.density) .- g.density) .* g.weights
  g.density ./= sum(g.density)
  g
end
function Adaptive(p::Int, F::Function, ::Type{T}, n::Int = 10_000, q::QuadratureRule = GenzKeister, l::Int = 6) where {F <: Function, T}
  ab = Adaptive(Dict{SVector{p, Int},T}(), AdaptiveBuildInit(p, q, F, l)...)
  construct!(ab, n)
  ab
end
function build(::Type{<: Adaptive{q, p, F, T}}, p::Int, F::Function, ::Type{T}, n::Int = 10_000, q::QuadratureRule = GenzKeister, l::Int = 6) where {}

end
gts(::Type{<:tsd{T}}) where {T} = tsd{T,Float64}(zero(T), [1.2, 5.3, 94.2])


function get!(GV::GridVessel{p, q, B}, i::Tuple{Int, Int}, F::Function where {p, q <: QuadratureRule}) where {B <: GridBuild}
  if haskey(GV.grids, i)
    return eval_grid(GV.grids[i], F)
  else
    grid_store, grid_out = build()
    GV.grids[i] = grid_store
    grid_out
  end
end
function get!(GV::GridVessel{p, q, B <: GridBuild}, i::Tuple{Int, Int}, F::Function, ::Type{T} where {p, q <: QuadratureRule}) where {B <: GridBuild, T}
  if haskey(GV.grids, i)
    return eval_grid_cache(GV.grids[i], F, T)
  else
    grid_store, grid_out = build(B)
    GV.grids[i] = grid_store
    grid_out
  end
end
function get!(GV::GridVessel{p, q, B <: aPrioriBuild}, i::Tuple{Int, Int},)
  if haskey(GV.grids, i)
    B()
  else
    something different.
  end
end
function get!(Grid::GridVessel{p, q, B <: AdaptiveBuild})
end

function FlattenedGrid(p::Int, l::Int, ::Type{q}, seq::Vector{Int} = default(q)) where {q}
  grid = NestedGrid(p, q, seq)
  smolyak!(grid, l)
  FlattenedGrid(grid)
end

function FlattenedGrid(Grid::NestedGrid{p,q}) where {p,q<:NestedQuadratureRule}
  n = length(Grid.grid)
  nodes = Array{Float64,2}(p,n)
  weights = Vector{Float64}(n)
  for (i, j) ∈ enumerate(keys(Grid.grid))
    weights[i] = Grid.grid[j]
    nodes[:,i] .= get_node.(Grid, j)
  end
  FlattenedGrid(nodes, weights, q)
end
FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, ::Type{GenzKeister}) = FlattenedGrid(n, w, - vec(sum(n .^ 2, [1])), GenzKeister)
FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, ::Type{KronrodPatterson}) = FlattenedGrid(n, w, zeros(w), KronrodPatterson)
function FlattenedGrid(n::Array{Float64,2}, w::Vector{Float64}, b::Vector{Float64}, ::Type{q}) where {q}
  FlattenedGrid{q}(n, w, b, similar(b))
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
