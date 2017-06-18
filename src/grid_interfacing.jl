struct SplitWeights{q<: QuadratureRule}
  nodes_positive::Array{Float64,2}
  nodes_negative::Array{Float64,2}
  weights_positive::Array{Float64,1}
  weights_negative::Array{Float64,1}
  baseline_weights_positive::Array{Float64,1}
  baseline_weights_negative::Array{Float64,1}
  weight_sum::Array{Float64,1}
  total_count::Int64
end

function SplitWeights(p::Int64, total_count::Int64, negative_count::Int64)
  positive_count = total_count - negative_count
  SplitWeights(Array{Float64,2}(p, positive_count), Array{Float64,2}(p, negative_count), Array{Float64,1}(positive_count), Array{Float64,1}(negative_count), Array{Float64,1}(positive_count), Array{Float64,1}(negative_count), Array{Float64,1}(2), total_count)
end

function SplitWeights{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q})
  nodes = maximum(values(Grid.Qs)).n
#  complete_grid = sum(values(Grid.Δ_prods))
  G = SplitWeights{q}(p, length(complete_grid), sum([i < 0 for i ∈ values(complete_grid)]))
  nc = 0; pc = 0
  for i ∈ keys(Grid.grid)
    w = complete_grid[i]
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

function set_baseline!(Grid::SplitWeights{GenzKeister})
  Grid.baseline_weights_positive .= - sum(Grid.nodes_positive .^ 2, [1])[1,:]
  Grid.baseline_weights_negative .= - sum(Grid.nodes_negative .^ 2, [1])[1,:]
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end
function set_baseline!(Grid::SplitWeights{KronrodPatterson})
  Grid.weight_sum[1] = sum(Grid.weights_positive)
  Grid.weight_sum[2] = sum(Grid.weights_negative)
end
