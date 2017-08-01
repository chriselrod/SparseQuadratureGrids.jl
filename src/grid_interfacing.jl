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

struct GridContainer{p,q}
  grids::Dict{Tuple{Int,Int},FlattenedGrid{q}}
  U::Array{Float64,2}
  l::Int
  seq::Vector{Int}
  mats::Dict{Int,Array{Float64,2}}
end

function getindex(g::GridContainer{p, q} where {p}, i::Int) where {q}
#  out = get!(g.grids, i, FlattenedGrid(i, g.l, q, g.seq))
#  copy!(out.density, out.baseline_weights)
#  out
  get!(g.grids, i, FlattenedGrid(i, g.l, q, g.seq))
end
function mats(g::GridContainer{p,q} where {q}, r::Int) where {p}
  get!(g.mats, r, Array{Float64,2}(p,r))
end
GridContainer(p::Int,l::Int = 6,::Type{q} = GenzKeister,seq::Vector{Int} = default(q)) where {q} = GridContainer{p,q}(Dict{Int,FlattenedGrid{q}}(), zeros(p,p), l, seq, Dict{Int,Array{Float64,2}}())

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
