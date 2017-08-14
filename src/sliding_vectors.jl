abstract type SlidingVector{T,P <: AbstractArray{T}} <: AbstractVector{T} end
mutable struct SlidingVec{T,P} <: SlidingVector{T, P}
    parent::P
    offset::Int
    l::Int
end
struct ChildVec{T,P} <: SlidingVector{T,P <: SlidingVector{T}}
    parent::P
    offset::Int
    l::Int
end
Base.IndexStyle(::Type{<: SlidingVector}) = IndexLinear()
function Base.getindex(V::SlidingVector, i::Int)
    @inbounds r = V.parent[V.offset + i]
    r
end
function Base.setindex!(V::SlidingVector, x, i::Int)
    @inbounds V.parent[V.offset + i] = x
    V
end

function Base.getindex(V::ChildVec, i::Int)
    @inbounds r = V.parent.parent[V.offset + V.parent.offset + i]
    r
end
function Base.setindex!(V::ChildVec, x, i::Int)
    @inbounds V.parent.parent[V.offset + V.parent.offset + i] = x
    V
end
Base.size(A::SlidingVector) = (A.l, )
Base.length(A::SlidingVector) = A.l
function SlidingVec(n::P, i::Int = 0, l::Int) where {T, P <: AbstractArray{T}}
    SlidingVec{T, P}(n, i, l)
end
function SlidingVec(n::P, i::Int = 0, l = size(n,2)) where {T, P <: AbstractArray{T,2}}
    SlidingVec{T, P}(n, i, l)
end
function SlidinvVec(n::P, i::UnitRange{Int}) where {T, P <: AbstractArray{T,2}}
    SlidingVec{T, P}(n, i[1], length(i))
end
function view(A::SlidingVec{T, P}, i::UnitRange{Int}) where {T, P}
    ChildVec{T, P}(A, i[1], length(i))
end
function slide!(A::SlidingVec)
    A.offset += A.l
end
function reset!(A::SlidingVec)
    A.offset = 0
end
