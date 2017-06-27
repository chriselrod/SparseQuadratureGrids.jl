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

function SplitWeights(p::Int64, total_count::Int64, negative_count::Int64, q::DataType)
  positive_count = total_count - negative_count
  SplitWeights{p,q}(Array{Float64,2}(p, positive_count), Array{Float64,2}(p, negative_count), Array{Float64,1}(positive_count), Array{Float64,1}(negative_count), Array{Float64,1}(positive_count), Array{Float64,1}(negative_count), Array{Float64,1}(2), total_count)
end

function SplitWeights{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q})
  nodes = maximum(values(Grid.Qs)).n
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

function SplitWeights{q <: QuadratureRule}(p::Int, l::Int, ::Type{q})
  grid = NestedGrid(p, q)
  smolyak!(grid, l)
  SplitWeights(grid)
end
function SplitWeights{q <: QuadratureRule}(p::Int, l::Int, ::Type{q}, seq::Vector{Int})
  grid = NestedGrid(p, q, seq)
  smolyak!(grid, l)
  SplitWeights(grid)
end

function set_baseline!{p}(Grid::SplitWeights{p, GenzKeister})
  Grid.baseline_weights_positive .= - sum(Grid.nodes_positive .^ 2, [1])[1,:]
  Grid.baseline_weights_negative .= - sum(Grid.nodes_negative .^ 2, [1])[1,:]
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end
function set_baseline!{p}(Grid::SplitWeights{p, KronrodPatterson})
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end
