abstract type SlidingVector{T,P <: AbstractArray{T}} <: AbstractVector{T} end
mutable struct SlidingVec{T,P} <: SlidingVector{T, P}
    parent::P
    offset::Int
    l::Int
end
struct ChildVec{T,P<: SlidingVector{T}} <: SlidingVector{T,P}
    parent::P
    offset::Int
    l::Int
end
struct SlidingVecFun{V, F}
    f::F
    v::V
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
function SlidingVec(n::P, i::Int = 0, l::Int = 1) where {T, P <: AbstractArray{T}}
    SlidingVec{T, P}(n, i, l)
end
function SlidingVec(n::P, i::Int = 0, l::Int = size(n,1)) where {T, P <: AbstractMatrix{T}}
    SlidingVec{T, P}(n, i, l)
end
function SlidingVec(n::P, i::UnitRange{Int}) where {T, P <: AbstractArray{T,2}}
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
function set!(A::SlidingVecFun, M::AbstractMatrix)
    A.v.parent = M
    A.v.offset = - A.v.l
end
function eval!(A::SlidingVecFun)
    A.v.offset += A.v.l
    A.f()
end
function (A::SlidingVecFun)(v::AbstractVecOrMat)
    A.v.parent[1+A.v.offset:length(v)+A.v.offset] .= v
    A.f()
end
function (A::SlidingVecFun)(v::SVector{p}) where p
    @inbounds for (i, v_i) ∈ enumerate(v)
        A.v.parent[i + A.v.offset] = v_i
    end
    A.v.parent[1+A.v.offset:A.v.l+A.v.offset] .= v
    A.f()
end
@inline (A::SlidingVecFun)() = eval!(A)
@inline eval!(A::SlidingVecFun, v::AbstractVector) = A(v)
#@inline (A::SlidingVecFun)(i::Int) = getindex(A, i)
function getindex(A::SlidingVecFun, i::Int)
    A.v.offset = A.v.l*(i-1)
    A.f()
end
function setindex!(A::SlidingVecFun, v::AbstractVector, i::Int)
    A.v.offset = A.v.l*(i-1)
    A.v.parent[1+A.v.offset:A.v.l+A.v.offset] .= v
    A.f()
end
function SlidingVecFun(f::Function, M::AbstractMatrix)
    l = size(M, 1)
    v = SlidingVec(M, -l, l)
    g() = f(v)
    SlidingVecFun(g, v)
end

struct randIter
    n::Int
end
function Base.start(A::SlidingVecFun)
    A.v.offset = 0
    A.f()
end
function Base.next(A::SlidingVecFun, state)
    eval!(A), state + 1
end
function Base.done(A::SlidingVecFun, state)
    state >= size(A.v.parent, 2)
end
function Base.eltype(A::SlidingVecFun)
    Core.Inference.return_type(A.f, ())
end
function Base.endof(A::SlidingVecFun)
    eval!(A, size(A.v.parent, 2))
end
Base.length(A::SlidingVecFun) = size(A.v.parent, 2)
Base.size(A::SlidingVecFun) = ( size(A.v.parent, 2), )
