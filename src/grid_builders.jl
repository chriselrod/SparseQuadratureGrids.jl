abstract type AdaptiveBuild{p, q, F} end
struct AdaptiveBuilder{p, q, F} <: AdaptiveBuild{p, q, F}
  F_cache::Dict{SVector{p,Int}, Float64}
  NeighborError::Dict{SVector{p,Int},Float64}
  Neighbors::Dict{SVector{p, Int}, Δprod{p}}
  Grid::NestedGrid{p, q}
  e::SVector{p, SVector{p,Int}}
  f::F
  l::Int #Max depth. 6 for GenzKeister & KronrodPatterson. Can be arbitrarilly high for GaussHermite.
end
struct AdaptiveBuilderCache{p, q, F, T} <: AdaptiveBuild{p, q, F}
  cache::Dict{SVector{p, Int}, T}
  F_cache::Dict{SVector{p,Int}, Float64}
  NeighborError::Dict{SVector{p,Int},Float64}
  Neighbors::Dict{SVector{p, Int}, Δprod{p}}
  Grid::NestedGrid{p, q}
  e::SVector{p, SVector{p,Int}}
  f::F
  l::Int
end

function AdaptiveBuildInit(p::Int, q::QuadratureRule, F::Function, l::Int = 6)
  Dict{SVector{p,Int}, Float64}(), Dict{SArray{p,Int},Float64}(), Dict{SVector{p, Int64}, Δprod{p}}(), NestedGrid(p, q), initialize_e(p), F, l
end

function AdaptiveBuilder(p::Int, q::QuadratureRule, F::Function, l::Int = 6)
  AdaptiveBuilder{p,q,F}(AdaptiveBuildInit(p, q, F, l)...)
end

function AdaptiveBuilder(p::Int, F::Function, n = 10_000, q::QuadratureRule = GenzKeister, l::Int = 6)
  ab = AdaptiveBuilder(p, q, F, l)
  construct!(ab, n)
  ab
end
function AdaptiveBuilderCache(p::Int, F::Function, ::Type{T}, n = 10_000, q::QuadratureRule = GenzKeister, l::Int = 6) where {T}
  ab = AdaptiveBuilderCache(Dict{SVector{p, Int},T}(), AdaptiveBuildInit(p, q, F, l)...)
  construct!(ab, n)
  ab
end

function e_j(p::Int, j::Int)
  out = zeros(Int, p)
  out[j] = 1
  SVector{p}(out)
end
function initialize_e(p::Int)
  SVector{p}(e_j.(p, 1:p))
end

function add_neighbors!(ab::AdaptiveBuild{p,q,F}, Δ_prod::SVector{p,Int}) where {p,q,F}
  for i ∈ 1:p
    Δ_prod_i = Δ_prod + ab.e[i]
    if haskey(ab.NeighborError, Δ_prod_i) || Δ_prod_i[i] > ab.l
      continue
    end
    neighbor_valid = true
    for j ∈ 1:p
      if ( Δprod[j] == 0 ) || ( i == j )
        continue
      elseif !haskey(ab.Grid.Δ_prods, Δprod_i - ab.e[j])
        neighbor_valid = false
        break
      end
    end
    if neighbor_valid
      add_neighbor!(ab, Δ_prod)
    end
  end
end
function add_neighbor!(ab::AdaptiveBuild{p,q,F}, Δ_prod::SVector{p,Int}) where {p,q,F}
  ab.Neighbors[Δ_prod] = calc_Δ_prod(Grid, Δ_prod)
  ab.NeighborError[Δ_prod] = Δ_prod_error!(ab, Δ_prod)
end



function construct!(ab::AdaptiveBuild{p,q,F}, n::Int) where {p,q,F}
  Δ_prod = @SVector zeros(Int, p)
  ab.Grid.Δ_prods[Δ_prod] = calc_Δ_prod!(ab.Grid, Δ_prod)
  while length(ab.Grid) < n
    add_neighbors!(ab, Δ_prod)
    Δ_prod = reduce((i,v) ->  i[2] > v[2] ? i : v, ab.NeighborError)[1]
    ab.Grid.Δ_prods[Δ_prod] = pop!(ab.Neighbors, Δ_prod)
    expand_grid!(ab.Grid, Δ_prod)
    delete!(ab.NeighborError, Δ_prod)
  end
end

function Δ_prod_error!(ab::AdaptiveBuild{p,q,F}, Δ_prod::SVector{p,Int}) where {p,q,F}
  out = 0.0
  for i ∈ keys(ab.Neighbors[Δ_prod].w)
    out += ab.Neighbors[Δ_prod].w[i] * eval_f!(ab, i)
  end
  out
end

function eval_f(ab::AdaptiveBuilder{p,q,F}, i::SVector{p,Int}) where {p,q,F}
  get!(() -> ab.f(get_node.(ab.Grid, i)), ab.F_cache, i)
end
function eval_f(ab::AdaptiveBuilderCache{p,q,F,T}, i::SVector{p,Int}) where {p,q,F,T}
  get!(() -> cache_f(ab, i), ab.F_cache, i)
end
function cache_f(ab::AdaptiveBuilderCache{p,q,F,T}, i::SVector{p,Int}) where {p,q,F,T}
  res, to_cache = ab.f(get_node.(ab.Grid, i))
  ab.cache[i] = to_cache
  res
end

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
