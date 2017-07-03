abstract type QuadratureRule end

abstract type NestedQuadratureRule <: QuadratureRule end
abstract type UnitNestedRule <: NestedQuadratureRule end
abstract type UnconstrainedNestedRule <: NestedQuadratureRule end

struct KronrodPatterson <: UnitNestedRule
#  nw::Dict{Int64,Tuple{Float64,Float64}}
  n::Dict{Int64,Float64}
  w::Dict{Int64,Float64}
  l::Int64
end
function KronrodPatterson(l::Int)
  dict = load(Pkg.dir("SparseQuadratureGrids")*"/src/rules/KronrodPatterson/l_"*string(l)*".jld")
  KronrodPatterson(dict["n"], dict["w"], l)
end
struct GenzKeister <: UnconstrainedNestedRule
#  nw::Dict{Int64,Tuple{Float64,Float64}}
  n::Dict{Int64,Float64}
  w::Dict{Int64,Float64}
  l::Int64
end
function GenzKeister(l::Int)
  dict = load(Pkg.dir("SparseQuadratureGrids")*"/src/rules/GenzKeister/l_"*string(l)*".jld")
  GenzKeister(dict["n"], dict["w"], l)
end
function delay_sequence(base_seq)
  seq = [1]
  pol_acc = 3
  seq_i = 2
  approx_pol_ac(x::Int) = round(Int, 3x/2+1/2)
  for i ∈ 2:length(base_seq)
    while approx_pol_ac(base_seq[i]) >= pol_acc
      push!(seq, base_seq[i])
      pol_acc += 2
    end
  end
  seq
end
function default(q::DataType)
  if q == GenzKeister
    return [1, 3, 9, 19, 35, 103]
    #return [1, 3, 3, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 19, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103]
  elseif q == KronrodPatterson
    return [1, 3, 7, 15, 31, 63]
    #return [1, 3, 3, 7, 7, 7, 15, 15, 15, 15, 15, 15, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63]
  else
    throw("Default unimplemented for grid type " * string(q) * ".")
  end
end

function Base.isless{T <: NestedQuadratureRule}(x::T, y::T)
  x.l < y.l
end

struct Δ{T <: NestedQuadratureRule}
  n::Dict{Int64,Float64}
  w::Dict{Int64,Float64}
end


function add_a_to_b{T, R <: Real}(a::Dict{T,R}, b::Dict{T,R})
  Dict{T,R}(i => b[i] + get(a, i, zero(R)) for i ∈ keys(b))
end
function Base.:-{T <: NestedQuadratureRule}(x::T, y::T)
  if x.l >= y.l
    return Δ{T}(x.n, Dict{Int64, Float64}(i => x.w[i] - get(y.w, i, 0.0) for i ∈ keys(x.w)))
  else
    return Δ{T}(y.n, Dict{Int64, Float64}(i => get(x.w,i,0.0) - y.w[i] for i ∈ keys(y.w)))
  end
end
function Base.:+{T, R <: Real}(x::Dict{T,R}, y::Dict{T,R})
    Dict{T, R}(i => get(x, i, zero(R)) + get(y, i, zero(R)) for i ∈ keys(x) ∪ keys(y))
end
function Base.:+{T, R <: Real}(x::Dict{T,R}...)
  Dict{T, R}(i => sum(get.(x, [i], zero(R))) for i ∈ union(keys.(x)...))
end
function Base.:+{T <: NestedQuadratureRule}(x::T, y::T)
  x.l > y.l ? Δ{T}(x.n, add_a_to_b(y.w, x.w)) : Δ{T}(y.n, add_a_to_b(x.w, y.w))
end
function Base.:+{T <: NestedQuadratureRule}(x::T...)
  max_q = maximum(x)
  Δ{T}(max_x.n, Dict{Int64, Float64}(i => sum(get.(x, [i], zero(R))) for i ∈ keys(max_x.w)))
end
