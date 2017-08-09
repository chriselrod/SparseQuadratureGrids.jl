abstract type GridBuild{q <: QuadratureRule} end
abstract type aPrioriBuild{q, p, F} <: GridBuild{q} end
abstract type AdaptiveBuild{q, p, F} <: GridBuild{q} end
struct AdaptiveRaw{q, p, F} <: AdaptiveBuild{q, p, F}
  F_cache::Dict{SVector{p,Int}, Float64}
  node_cache::Dict{SVector{p,Int}, Vector{Float64}}
  baseline_cache::Dict{SVector{p,Int}, Float64}
  NeighborError::Dict{SVector{p,Int}, Float64}
  Neighbors::Dict{SVector{p, Int}, Δprod{p}}
  Grid::NestedGrid{p, q}
  e::SVector{p, SVector{p,Int}}
  f::F
  l::Int #Max depth. 6 for GenzKeister & KronrodPatterson. Can be arbitrarilly high for GaussHermite.
end
struct Adaptive{q, p, F, T} <: AdaptiveBuild{q, p, F}
  cache::Dict{SVector{p, Int}, T}
  F_cache::Dict{SVector{p,Int}, Float64}
  node_cache::Dict{SVector{p,Int}, Vector{Float64}}
  baseline_cache::Dict{SVector{p,Int}, Float64}
  NeighborError::Dict{SVector{p,Int}, Float64}
  Neighbors::Dict{SVector{p, Int}, Δprod{p}}
  Grid::NestedGrid{p, q}
  e::SVector{p, SVector{p,Int}}
  f::F
  l::Int
end
struct SmolyakRaw{q, p, F} <: aPrioriBuild{q, p, F} end
struct Smolyak{q, p, F, T} <: aPrioriBuild{q, p, F} end
RawBuild{q, p, F} = Union{SmolyakRaw{q, p, F}, AdaptiveRaw{q, p, F}}
SmolyakBuild{q, p, F} = Union{SmolyakRaw{q, p, F}, Smolyak{q, p, F, T} where T}

function e_j(::Type{Val{p}}, j::Int) where p
  out = zeros(Int, p)
  out[j] = 1
  SVector{p}(out)
end
function initialize_e(::Type{Val{p}}) where p
  SVector{p}(e_j.(Val{p}, 1:p))
end

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
  Δ_prod = @SVector ones(Int, p)
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
  out
end
function get_node(ab::AdaptiveBuild{q,p,F} where {q,F}, i::SVector{p,Int}) where p
  get!(() -> get_node.(ab.Grid, i), ab.node_cache, i)
end

function eval_f!(ab::AdaptiveRaw{GenzKeister,p,F} where F, i::SVector{p,Int}) where p
  get!(() -> cache_f(ab, i), ab.F_cache, i)
end
function eval_f!(ab::AdaptiveRaw{KronrodPatterson,p,F} where F, i::SVector{p,Int}) where p
  get!(() -> exp(ab.f(get_node(ab, i))), ab.F_cache, i)
end
function eval_f!(ab::Adaptive{q,p,F,T} where {q,F,T}, i::SVector{p,Int}) where p
  get!(() -> cache_f(ab, i), ab.F_cache, i)
end
function cache_f(ab::AdaptiveRaw{GenzKeister,p,F} where F, i::SVector{p,Int}) where p
  node = get_node(ab, i)
  ab.baseline_cache[i] = sum(node .^ 2)
  exp(ab.f(node) + ab.baseline_cache[i])
end
function cache_f(ab::Adaptive{GenzKeister,p,F,T} where {F,T}, i::SVector{p,Int}) where p
  node = get_node(ab, i)
  ab.baseline_cache[i] = sum(node .^ 2)
  res, ab.cache[i] = ab.f(node)
  exp(res + ab.baseline_cache[i])
end
function cache_f(ab::Adaptive{KronrodPatterson,p,F,T} where {F,T}, i::SVector{p,Int}) where p
  res, ab.cache[i] = ab.f(get_node(ab, i))
  exp(res)
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
