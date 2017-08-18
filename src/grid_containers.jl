#Why the p in grid interfaces?
abstract type GridInterface{q <: QuadratureRule} end
abstract type GridVessel{q, B <: GridBuild, K, p}  <: GridInterface{q} end
abstract type CollapsedGrid{q <: QuadratureRule} end
struct FlatGrid{q,N,C} <: CollapsedGrid{q}
    nodes::N
    cache::C
    weights::Vector{Float64}
    density::Vector{Float64}
    baseline::Vector{Float64}
end
struct MatrixCache{p}
  d::Dict{Int,Array{Float64,2}}
end
struct DynamicGridVessel{q, B, K, p, N, C} <: GridVessel{q, B, K, p}
    grids::Dict{K, FlatGrid{q,N,C}}
    U::Matrix{Float64}
    mats::MatrixCache{p}
end
struct MatrixVecSVec{p, T <: Real}
    M::Matrix{T}
    S::Vector{SVector{p,T}}
end
