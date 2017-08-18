abstract type GridBuild{q <: QuadratureRule} end
abstract type aPrioriBuild{q, p, F} <: GridBuild{q} end
abstract type AdaptiveBuild{q, p, F, A, B} <: GridBuild{q} end
struct AdaptiveRaw{q, p, F, A, B, d} <: AdaptiveBuild{q, p, F, A, B}
  F_cache::Dict{SVector{p,Int}, Float64}
  cache::Dict{SVector{p,Int}, SVector{d,Float64}}
  node_cache::Dict{SVector{p,Int}, SVector{p,Float64}}
  baseline_cache::Dict{SVector{p,Int}, Float64}
  NeighborError::Dict{SVector{p,Int}, Float64}
  Neighbors::Dict{SVector{p, Int}, Δprod{p}}
  Grid::NestedGrid{p, q}
  e::SVector{p, SVector{p,Int}}
  f::F
  l::Int
  μ::A
  U::B
end
struct Adaptive{q, p, F, A, B, T} <: AdaptiveBuild{q, p, F, A, B}
  cache::Dict{SVector{p, Int}, T}
  F_cache::Dict{SVector{p,Int}, Float64}
  node_cache::Dict{SVector{p,Int}, SVector{p,Float64}}
  baseline_cache::Dict{SVector{p,Int}, Float64}
  NeighborError::Dict{SVector{p,Int}, Float64}
  Neighbors::Dict{SVector{p, Int}, Δprod{p}}
  Grid::NestedGrid{p, q}
  e::SVector{p, SVector{p,Int}}
  f::F
  l::Int
  μ::A
  U::B
end
struct SmolyakRaw{q, p, F} <: aPrioriBuild{q, p, F} end
struct Smolyak{q, p, F, T} <: aPrioriBuild{q, p, F} end
RawBuild{q, p, F} = Union{SmolyakRaw{q, p, F}, AdaptiveRaw{q, p, F}}
CacheBuild{q, p, F} = Union{Adaptive{q, p, F}, Smolyak{q, p, F}}
SmolyakBuild{q, p, F} = Union{SmolyakRaw{q, p, F}, Smolyak{q, p, F, T} where T}
default(::Type{AB}) where {AB <: AdaptiveBuild} = 1_000
default(::Type{AP}) where {q, AP <: aPrioriBuild{q}} = default(q)

function Base.sizehint!(ab::AdaptiveBuild, n::Int)
    sizehint!(ab.cache, n)
    sizehint!(ab.F_cache, n)
    sizehint!(ab.node_cache, n)
    sizehint!(ab.baseline_cache, n)
    sizehint!(ab.Grid.grid.w, n)
end

initialize_e(::Type{Val{p}}) where p = SVector{p}(svectors(eye(Int, p), Val{p}))

function add_neighbors!(ab::AdaptiveBuild{q,p,F} where {q,F}, Δ_prod::SVector{p,Int}) where p
  for i ∈ 1:p
    Δ_prod_i = Δ_prod + ab.e[i]
    if haskey(ab.NeighborError, Δ_prod_i) || Δ_prod_i[i] > ab.l
      continue
    end
    neighbor_valid = true
    for j ∈ 1:p
      if ( Δ_prod[j] == 1 ) || ( i == j )
        continue
      elseif !haskey(ab.Grid.Δ_prods, Δ_prod_i - ab.e[j])
        neighbor_valid = false
        break
      end
    end
    if neighbor_valid
      add_neighbor!(ab, Δ_prod_i)
    end
  end
end
function add_neighbor!(ab::AdaptiveBuild{q,p,F} where {q,F}, Δ_prod::SVector{p,Int}) where p
  ab.Neighbors[Δ_prod] = calc_Δ_prod(ab.Grid, Δ_prod)
  ab.NeighborError[Δ_prod] = Δ_prod_error!(ab, Δ_prod)
end



function construct!(ab::AdaptiveBuild{q,p,F} where {q,F}, n::Int) where p
    sizehint!(ab, n)
    Δ_prod = SVector{p}(ones(Int, p))
    ab.Grid.Δ_prods[Δ_prod] = calc_Δ_prod!(ab.Grid, Δ_prod)
    while length(ab.Grid) < n
        add_neighbors!(ab, Δ_prod)
        Δ_prod = reduce((i,v) ->  i[2] > v[2] ? i : v, ab.NeighborError)[1]
        ab.Grid.Δ_prods[Δ_prod] = pop!(ab.Neighbors, Δ_prod)
        expand_grid!(ab.Grid, Δ_prod)
        delete!(ab.NeighborError, Δ_prod)
    end
end

function Δ_prod_error!(ab::AdaptiveBuild{q,p,F} where {q,F}, Δ_prod::SVector{p,Int}) where p
  out = 0.0
  for (i, w) ∈ ab.Neighbors[Δ_prod].w
    out += w * eval_f!(ab, i)
  end
  abs(out)
end
function get_node(ab::AdaptiveBuild{GenzKeister,p,F} where {F}, i::SVector{p,Int}) where p
  get!(() -> get_node.(ab.Grid, i), ab.node_cache, i)
end
function get_node(ab::AdaptiveBuild{KronrodPatterson,p,F} where {F}, i::SVector{p,Int}) where p
    get!(() -> calc_node(ab, i), ab.node_cache, i)
end
function calc_node(ab::AdaptiveBuild{KronrodPatterson}, i::SVector{p,Int}) where p
    snode = get_node.(ab.Grid, i)
    if !haskey(ab.baseline_cache, i)
        ab.baseline_cache[i] = sum(log_sigmoid_jacobian, snode)
    end
    sigmoid.(snode)
end
function transformed_node(ab::AdaptiveRaw, i::SVector{p,Int}, v::SVector{p,Float64}) where p
    get!(() -> ab.μ + ab.U * v, ab.cache, i)
end

function eval_f!(ab::AdaptiveBuild, i::SVector)
  get!(() -> cache_f(ab, i), ab.F_cache, i)
end
function cache_f(ab::AdaptiveRaw, i::SVector)
  node = get_node(ab, i)
  exp( ab.f(transformed_node(ab, i, node)) + get_baseline!(ab, i, node) )
end
function cache_f(ab::Adaptive, i::SVector)
  node = get_node(ab, i)
  res, ab.cache[i] = ab.f(ab.μ + ab.U * node)
  exp( res + get_baseline!(ab, i, node) )
end
@inline function get_baseline!(ab::AdaptiveBuild{KronrodPatterson}, i::SVector, ::AbstractVector)
    get!( () -> sum(log_sigmoid_jacobian, get_node.(ab.Grid, i)), ab.baseline_cache, i )
end
@inline function get_baseline!(ab::AdaptiveBuild{GenzKeister}, i::SVector, x::AbstractVector)
    get!( () -> sum(abs2, x), ab.baseline_cache, i )
end
@inline sigmoid(x::Real) = log( (1 + x)/(1 - x) )
@inline log_sigmoid_jacobian(x::Real) = - log( 1 - x^2 )

@generated function smolyak!(Grid::NestedGrid{p,q}, l::Int) where {p,q<:NestedQuadratureRule}
  quote
    eval(parse("j_"*string($p+1)*" = "*string(l)))
    eval(parse("s_"*string($p+1)*" = 0"))
    @nloops $p i d -> begin
      1:j_{d+1}
    end d -> begin
      s_d = s_{d+1} + i_d - 1
      j_d = l - s_d
    end begin
      Δ = SVector{p}((@ntuple $p i))
      Grid.Δ_prods[Δ] = calc_Δ_prod!(Grid, Δ)
    end
  end
end

@generated function l2!(Grid::NestedGrid{p,q}, l::Int) where {p,q<:NestedQuadratureRule}
  quote
    eval(parse("j_"*string($p+1)*" = "*string(l)))
    eval(parse("s_"*string($p+1)*" = 0"))
    r2 = (l-1)^2
    @nloops $p i d -> begin
      1:j_{d+1}
    end d -> begin
      s_d = s_{d+1} + (i_d - 1)^2
      j_d = floor(Int, √(r2 - s_d) )
    end begin
      Δ = SVector{p}((@ntuple $p i))
      Grid.Δ_prods[Δ] = calc_Δ_prod!(Grid, Δ)
    end
  end
end


@generated function lx!(Grid::NestedGrid{p,q}, l::Int, x::Real) where {p,q<:NestedQuadratureRule}
  quote
    eval(parse("j_"*string($p+1)*" = "*string(l)))
    eval(parse("s_"*string($p+1)*" = 0"))
    rx = (l-1)^x
    inv_x = 1/x
    @nloops $p i d -> begin
      1:j_{d+1}
    end d -> begin
      s_d = s_{d+1} + (i_d - 1)^x
      j_d = floor(Int, (rx - s_d)^inv_x )
    end begin
      Δ = SVector{p}((@ntuple $p i))
      Grid.Δ_prods[Δ] = calc_Δ_prod!(Grid, Δ)
    end
  end
end
