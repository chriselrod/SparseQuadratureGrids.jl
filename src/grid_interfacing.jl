function DynamicGridVessel(::Type{q},::Type{B},::Type{K},::Type{Val{p}}, ::Type{N}, ::Type{C}) where {q <: NestedQuadratureRule, B <: GridBuild, K, p, N, C}
  DynamicGridVessel{q,B,K,p,N,C}(Dict{K,FlatGrid{q,N,C}}(),Array{Float64,2}(p,p), MatrixCache(Val{p}))
end


MatrixCache(::Type{Val{p}}) where p = MatrixCache{p}(Dict{Int,Array{Float64,2}}())
MatrixCache(p::Int) = MatrixCache{p}(Dict{Int,Array{Float64,2}}())

function MatrixVecSVec(x::Matrix, ::Type{Val{p}}) where p
    MatrixVecSVec(x, svectors(x, Val{p}))
end
function MatrixVecSVec(x::Vector{SVector{p,T}} where T<:Real) where p
    MatrixVecSVec(mat(x), x)
end
function svectors(x::Matrix{T}, ::Type{Val{p}}) where {p,T}
    reinterpret(SVector{p,T}, x, (size(x,2),))
end
function mat(x::Vector{SVector{p,Float64}}) where p
    reinterpret(Float64, x, ( p, length(x) ) )
end
@inline convert(::Type{Vector{SVector{p,T}}}, x::Matrix{T}) where {p,T} = svectors(x, Val{p})
@inline convert(::Type{Matrix{T}}, x::Vector{SVector{p,T}} where p) where T = mat(x)
@inline convert(::Type{MatrixVecSVec{p, T}}, x::Matrix{T} ) where {p,T} = MatrixVecSVec(x, Val{p})
@inline convert(::Type{MatrixVecSVec{p, T}}, x::Vector{SVector{p,T}} ) where {p,T} = MatrixVecSVec(x)


### These functions update the cached grid in place.
#LogDensities' Model() function uses heuristics to choose Benchmark-based heurstics (ie, >13 ? dimensions) to choose whether iterating over SVectors or a BLAS call are more performant.
@inline update_cache_grid!(nodes::Matrix, cache::MatrixVecSVec, μ, U) = update_cache_grid(nodes, cache.M, μ, U)
function update_cache_grid!(nodes::Matrix, cache::Matrix, μ, U)
    cache .= μ
    Base.LinAlg.BLAS.gemm!('N','N', 1.0, U, nodes, 1.0, cache)
end
@inline update_cache_grid!(nodes::Vector{SVector{p,T}} where {p,T}, cache::MatrixVecSVec, μ, U) = update_cache_grid!( nodes, cache.S, μ, U )
function update_cache_grid!(nodes::Vector{SVector{p,T}}, cache::Vector{SVector{d,T}}, μ, U) where {p, d, T}
    U_sm = SMatrix{d,p}(U)
    μ_sv = SVector{d}(μ)
    @inbounds for (i, j) in enumerate(nodes)
        cache[i] = μ_sv + U_sm * j
    end
end


#A sliding view evaluates all the grid points.
#This saves on any reading/writing one may have to do to memory.
#The advantage may be negligble for functions of simple vectors, but could save a lot of memory rewriting / object construction for other data types, like those offered in ConstrainedParameters.
function eval_f!(cache::Matrix, weights::Vector, density::Vector, baseline::Vector, svf::SlidingVecFun)
    set!(svf, cache)
    @inbounds for i ∈ eachindex(density)
        density[i] = weights[i] * exp( baseline[i] + eval!(svf) )
    end
    density ./= sum(density)
end
@inline eval_f!(cache::MatrixVecSVec, weights::Vector, density::Vector, baseline::Vector, svf::SlidingVecFun) = eval_f!(cache.M, weights, density, baseline, svf)
@inline eval_f!(FG::FlatGrid, svf::SlidingVecFun) = eval_f!(FG.cache, FG.weights, FG.density, FG.baseline, svf)
function eval_grid!(FG::FlatGrid, svf::SlidingVecFun, μ::AbstractVector, U::AbstractMatrix)
    update_cache_grid!(FG.nodes, FG.cache, μ, U)
    eval_f!(FG, svf)
    FG
end
function eval_grid!(FG::FlatGrid, f::Function, μ::AbstractVector, U::AbstractMatrix, ::Type{T}) where T
    update_cache_grid!(FG.nodes, FG.cache, μ, U)
    eval_f!(FG, FG.cache, f, T)
end

@inline eval_f!(g::FlatGrid, cache::MatrixVecSVec{p,Float64}, f::Function, ::Type{T}) where {p,T} = eval_f!(g, cache.S, f, T)
function eval_f!(g::FlatGrid, cache::Vector{SVector{p,Float64}}, f::Function, ::Type{T}) where {p,T}
    param_cache = map( enumerate( cache ) ) do k
        i, j = k
        f_val, p_cache = f( j )
        g.density[i] = g.weights[i] * exp( g.baseline[i] + f_val )
        p_cache
    end
#    Vector{T}( length( g.weights ) )
#    @inbounds for (i,j) ∈ enumerate( cache )
#        f_val, param_cache[i] = f( j )
#        g.density[i] = g.weights[i] * exp( g.baseline[i] + f_val )
#    end
    g.density ./= sum(g.density)
    param_cache, g.density
end
function Adapt(::Type{Adaptive{q, p, F, T}}, f::F, l::Int, n::Int, μ::A, U::B) where {q,p, F <: Function, A, B, T}
    ab = Adaptive{q, p, F, A, B, T}(Dict{SVector{p, Int},T}(), Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int}, SVector{p, Float64}}(), Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int},Float64}(), Dict{SVector{p, Int64}, Δprod{p}}(), NestedGrid(Val{p}, q), initialize_e(Val{p}), f, l, μ, U)
    construct!(ab, n)
    ab
end
function Adapt(::Type{AdaptiveRaw{q, p, SVF}}, svf::SVF, l::Int, n::Int, μ::A, U::B) where {q,p, SVF <: SlidingVecFun, d, A <: SVector{d}, B}
    ab = AdaptiveRaw{q, p, SVF, A, B, d}(Dict{SVector{p,Int}, Float64}(),Dict{SVector{p,Int}, SVector{d,Float64}}(), Dict{SVector{p,Int}, SVector{p, Float64}}(), Dict{SVector{p,Int}, Float64}(), Dict{SVector{p,Int},Float64}(), Dict{SVector{p, Int64}, Δprod{p}}(), NestedGrid(Val{p}, q), initialize_e(Val{p}), svf, l, μ, U)
    construct!(ab, n)
    ab
end


Adaptive(::Type{Val{p}}, ::F, ::Type{T}, n::Int = 1_000, ::Type{q} = GenzKeister, l::Int = 6) where {q,p,F,T} = Adaptive(Adaptive{q, p, F, T}, l, n)

function density(ab::Adaptive{q, d, F, A, B, T} where {q, d, F, A, B}, ::Type{N} where N, ::Type{C} where C, ::Type{Val{p}} where p) where T
  n = length(ab.Grid.grid)
  density = Vector{Float64}(n)
  params = Vector{T}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    density[i] = ab.F_cache[j] * k
    params[i] = ab.cache[j]
  end
  density ./= sum(density)
  params, density
end
function density(ab::AdaptiveRaw{q, d, F} where F, ::Type{N}, ::Type{C}, ::Type{Val{p}}) where {q,p,d,N,C}
  n = length(ab.Grid.grid)
  nodes = Array{SVector{d,Float64}}(n)
  cache = Array{SVector{p,Float64}}(n)
  weights = Vector{Float64}(n)
  density = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    density[i] = ab.F_cache[j] * k
    nodes[i] = get_node(ab, j)
    cache[i] = ab.cache[j]
  end
  density ./= sum(density)
  FlatGrid{q,N,C}(convert(N,nodes), convert(C,cache), weights, density, Vector{Float64}(n))
end

 @noinline function build!(GV::DynamicGridVessel{q, B, K, p, N, C}, ::Type{Bc}, f, i::K, μ::Vector, U::Matrix, l::Int = 6) where {q, p, d, B <: AdaptiveBuild{q}, Bc <: AdaptiveBuild{q, d}, K, N, C}
  ab = Adapt(Bc, f, l, i[end], SVector{p}(μ), SizedArray{Tuple{p,d}}(U))
  GV.grids[i] = FlatGrid(ab, N, C, Val{p})
  for v ∈ values(ab.Neighbors)
    expand_grid!(ab.Grid, v)
  end
  density(ab, N, C, Val{p})
end
function build!(GV::DynamicGridVessel{q, B, K, p, N, C}, ::Type{Bc}, svf::SVF, i::K, μ::AbstractVector, U::AbstractMatrix) where {q, p, SVF <: SlidingVecFun, B <: aPrioriBuild{q}, Bc <: aPrioriBuild{q,p,SVF}, K, N, C}
    GV.grids[i] = grid = FlatGrid(Bc, ind(i), Val{p}, N, C)
    eval_grid!(grid, svf, μ, U)
end
function build!(GV::DynamicGridVessel{q, B, K, p, N, C}, ::Type{Bc}, f::F, i::K, μ::AbstractVector, U::AbstractMatrix) where {q, p, T, F<:Function, B <: aPrioriBuild{q}, Bc <: Smolyak{q,p,F,T}, K, N, C}
    GV.grids[i] = grid = FlatGrid(Bc, ind(i), Val{p}, N, C)
    eval_grid!(grid, f, μ, U, T)
end
@inline ind(i::Vector{Int}) = i
@inline ind(i::Tuple{Int,Vector{Int}}) = i[2]




#Oops, the ps for dispatch here should be pulled from R when R <: Static Rank.
#This would also allow for soemthing cleaner than "Dynamic Inds."
#Instead: given Val{d} ? use d : use i[1]
@inline function calc_grid!(Grid::DynamicGridVessel{q, B, K, p} where {q <: QuadratureRule}, i::K, svf::SVF, μ::Vector, U::Matrix) where {p, B <: RawBuild, K, SVF <: SlidingVecFun}
  build!(Grid, B{p, SVF}, svf, i, μ, U)
end
@inline function calc_grid!(Grid::DynamicGridVessel{q, B, K, p} where {q <: QuadratureRule}, i::K, svf::F, μ::Vector, U::Matrix, ::Type{T}) where {p, B <: CacheBuild, K, F <: Function, T}
  build!(Grid, B{p, F, T}, svf, i, μ, U)
end


#This is probably really bad, but these indices "K" are used by JointPosteriors when a dynamic dispatch on grid construction is necessary.
#In this case, these sets of tuples are generated for Dynamic Rank models, where the dimensionality is reduced by variable amounts via LDR.
#If one has an idea of how many parameters they wish to keep, "Fixed{p}" is better LDR because it allows us to avoid this dynamic dispatch.
DynamicInds = Union{Tuple{Int,Int,Int},Tuple{Int,Vector{Int}}}
function calc_grid!(Grid::DynamicGridVessel{q, B, K, p, N, C} where {q <: QuadratureRule}, i::K, svf::SVF, μ::Vector, U::Matrix) where {p, B <: RawBuild, K <: DynamicInds, SVF <: SlidingVecFun, N, C}
  build!(Grid, B{i[1], SVF}, svf, i, μ, U)::FlatGrid{q,N,C}
end
function calc_grid!(Grid::DynamicGridVessel{q, B, K, p, N, C} where {q <: QuadratureRule}, i::K, f::F, μ::Vector, U::Matrix, ::Type{T}) where {p, B <: CacheBuild, K <: DynamicInds, N, F, C, T}
  build!(Grid, B{i[1], F, T}, f, i, μ, U)::FlatGrid{q,N,C}
end

function eval_grid!( Grid::G, svf::SlidingVecFun, μ::Vector, U::Matrix, i::K ) where {q, B, K, G <: GridVessel{q, B, K}}
    haskey(Grid.grids, i) ? eval_grid!(Grid.grids[i], svf, μ, U) : calc_grid!(Grid, i, svf, μ, U)
end
function eval_grid!( Grid::G, f::Function, μ::Vector, U::Matrix, i::K, ::Type{MP} ) where {q, B, K, G <: GridVessel{q, B, K}, MP}
    haskey(Grid.grids, i) ? eval_grid!(Grid.grids[i], f, μ, U, MP) : calc_grid!(Grid, i, f, μ, U, MP)
end





@inline calc_baseline(::Type{GenzKeister}, n::Vector{SVector{p,T}}) where {p,T} = sum.(abs2, n)
@inline calc_baseline(::Type{KronrodPatterson}, n::Vector{SVector{p,T}}) where {p,T} = sum.(log_sigmoid_jacobian, n)

@inline transform_nodes(::Type{GenzKeister}, ::Type{N} where N, n) = n
@inline transform_nodes(::Type{KronrodPatterson}, ::Type{Matrix{Float64}}, n) = sigmoid.(mat(n))
function transform_nodes(::Type{KronrodPatterson}, ::Type{Vector{SVector{p,Float64}}}, n::Vector{SVector{p,Float64}}) where p
    @inbounds for (i, j) ∈ enumerate(n)
        n[i] = sigmoid.(j)
    end
    n
end
function FlatGrid(::Type{<:SmolyakBuild{q,d}}, seq::Vector{Int}, ::Type{p}, ::Type{N}, ::Type{C}) where {q,d,p,N,C}
  grid = NestedGrid(p, q, seq)
  smolyak!(grid, length(seq))
  FlatGrid(grid, p, N, C)
end
function FlatGrid(Grid::NestedGrid{d,q}, ::Type{p}, ::Type{N}, ::Type{C}) where {q,d,p,N,C}
  n = length(Grid.grid)
  nodes = Vector{SVector{d,Float64}}(n)
  weights = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(Grid.grid.w)
    weights[i] = k
    nodes[i] = get_node.(Grid, j)
  end
  base = calc_baseline(q, nodes)
  FlatGrid(transform_nodes(q, N, nodes), weights, base, Vector{Float64}(n), q, N, C, p, n)
end
function FlatGrid(ab::Adaptive{q,d,F}, ::Type{N}, ::Type{C}, ::Type{Val{p}}) where {q,d,p,F,N,C}
  n = length(ab.Grid.grid)
  nodes = Vector{SVector{d,Float64}}(n)
  weights = Vector{Float64}(n)
  baseline = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    nodes[i] = get_node(ab, j)
    baseline[i] = ab.baseline_cache[j]
  end
  FlatGrid(nodes, weights, baseline, Vector{Float64}(n), q, N, C, Val{p}, n)
end
function FlatGrid(ab::AdaptiveRaw{q,d,F}, ::Type{N}, ::Type{C}, ::Type{Val{p}}) where {q,d,p,F,N,C}
  n = length(ab.Grid.grid)
  nodes = Vector{SVector{d,Float64}}(n)
  weights = Vector{Float64}(n)
  baseline = Vector{Float64}(n)
  @inbounds for (i, (j,k)) ∈ enumerate(ab.Grid.grid.w)
    weights[i] = k
    nodes[i] = get_node(ab, j)
    baseline[i] = ab.baseline_cache[j]
  end
  FlatGrid(nodes, weights, baseline, Vector{Float64}(n), q, N, C, Val{p}, n)
end
@inline function FlatGrid(nodes::Vector{SVector{d,Float64}} where d, weights::Vector, baseline::Vector, density::Vector, ::Type{q}, ::Type{Matrix{Float64}}, ::Type{Matrix{Float64}}, ::Type{Val{p}}, n::Int) where {q, p}
    FlatGrid{q, Matrix{Float64}, Matrix{Float64}}(mat(nodes), Array{Float64,2}( p , n ), weights, density, baseline)
end
@inline function FlatGrid(nodes::N, weights::Vector, baseline::Vector, density::Vector, ::Type{q}, ::Type{N}, ::Type{C}, ::Type{Val{p}}, n::Int) where {q, p, N, C}
    FlatGrid{q, N, C}(nodes, convert(C, Vector{SVector{p,Float64}}(n)), weights, density, baseline)
end
@inline function FlatGrid(nodes::Vector{SVector{d,Float64}} where d, weights::Vector, baseline::Vector, density::Vector, ::Type{q}, ::Type{Matrix{Float64}}, ::Type{MatrixVecSVec{p,Float64}}, ::Type{Val{p}}, n::Int) where {q, p}
    FlatGrid{q, Matrix{Float64}, MatrixVecSVec{p,Float64}}(mat(nodes),
    MatrixVecSVec(Matrix{Float64}(p,n), Val{$p}), weights, density, baseline)
end

mats(mc::MatrixCache{p}, r) where p = get!(() -> Array{Float64,2}(p,r), mc.d, r)
mats(mc::MatrixCache{p}) where p = get!(() -> Array{Float64,2}(p,p), mc.d, p)
@inline mats(g::GridVessel, r::Int) = mats(g.mats, r)
@inline mats(g::GridVessel) = mats(g.mats)
