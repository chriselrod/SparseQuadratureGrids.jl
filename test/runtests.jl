using SparseQuadratureGrids
using Base.Test

# write your own tests here


Base.var(::Type{GenzKeister}) = 1/2
Base.var(::Type{KronrodPatterson}) = 1/3

function testSplitWeights(g::SplitWeights{p,q}, ec::Int) where {p,q}
  @testset begin
    @test sum(g.weight_sum) ≈ 1.0
    @test g.total_count == ec
    @test minimum(g.nodes_positive) < minimum(g.nodes_negative)
    @test maximum(g.nodes_positive) > maximum(g.nodes_negative)
    @test all(isapprox.(g.nodes_positive * g.weights_positive, 0.0, atol = 1e-14))
    @test all(isapprox.(g.nodes_negative * g.weights_negative, 0.0, atol = 1e-14))
    @test all(isapprox.(g.nodes_positive .^ 2 * g.weights_positive .+ g.nodes_negative .^ 2 * g.weights_negative, var(q)))
  end
end

dsgk = [1, 3, 3, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 19, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 103]
dskp = [1, 3, 3, 7, 7, 7, 15, 15, 15, 15, 15, 15, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 63]
testSplitWeights(SplitWeights(10,6,GenzKeister), 115645)
testSplitWeights(SplitWeights(4,length(dsgk),GenzKeister,dsgk), 253649)

testSplitWeights(SplitWeights(8,6,KronrodPatterson), 31745)
testSplitWeights(SplitWeights(5,length(dskp),KronrodPatterson,dskp), 1297727)
