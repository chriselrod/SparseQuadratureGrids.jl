struct Val{p} end
struct Δprod{p}
  w::Dict{SVector{p, Int64},Float64}
end
function Δprod(::Type{Val{p}}) where p
  Δprod{p}(Dict{SVector{p, Int64},Float64}())
end
function Base.:+(x::Δprod{p}, y::Δprod{p}) where p
  Δprod(x.w + y.w)
end
function Base.:+(x::Δprod{p}...) where p
  Δprod(sum(x.w))
end
function Base.setindex!(A::Δprod{p}, x::Float64, i::Int64...) where {p}
  A.w[SVector{p,Int64}(i)] = x
end
function Base.getindex(A::Δprod{p}, i::SVector{p,Int64}) where p
  A.w[i]
end
function Base.keys(A::Δprod)
  keys(A.w)
end
function Base.values(A::Δprod)
  values(A.w)
end
function Base.length(A::Δprod)
  length(A.w)
end


struct NestedGrid{p, q <: NestedQuadratureRule}
  Qs::Dict{Int64, q}
  Δs::Dict{Int64, Δ{q}}
  Δ_prods::Dict{SVector{p, Int64}, Δprod{p}}
  node_value::Dict{Int64,Float64}
  seq::Array{Int64,1}
  grid::Δprod{p}
#  w::Dict{SVector{p, Int64}, Float64}

end
function NestedGrid(::Type{Val{p}}, ::Type{q}) where {p,q<:NestedQuadratureRule}
  NestedGrid{p, q}(Dict{Int64,q}(), Dict{Int64,Δ{q}}(), Dict{SVector{p, Int64}, Δprod}(), Dict{Int64,Float64}(), default(q), Δprod(p))
end
function NestedGrid(::Type{Val{p}}, ::Type{q}, seq::Array{Int64,1}) where {p,q<:NestedQuadratureRule}
  NestedGrid{p, q}(Dict{Int64,q}(), Dict{Int64,Δ{q}}(), Dict{SVector{p, Int64}, Δprod}(), Dict{Int64,Float64}(), seq, Δprod(Val{p}))
end

#Type unstable calls.
NestedGrid(p::Int, ::Type{q}) = NestedGrid(Val{p}, q)
NestedGrid(p::Int, ::Type{q}, seq::Array{Int64,1}) = NestedGrid(Val{p}, q, seq)

function get_node(Grid::NestedGrid, i::Int)
  Grid.node_value[i]
end
function calc_Q!(Grid::NestedGrid{p, q}, i::Int) where {p, q <: NestedQuadratureRule}
  new_q = q(Grid.seq[i])
  if length(new_q.n) > length(Grid.node_value)
    merge!(Grid.node_value, new_q.n)
  end
  new_q
end
function get_Q!(Grid::NestedGrid{p, q}, i::Int) where {p, q <: NestedQuadratureRule}
  get!(() -> calc_Q!(Grid, i), Grid.Qs, i)
end
function calc_Δ!(Grid::NestedGrid{p,q}, i::Int) where {p, q <: NestedQuadratureRule}
  if i > 1
    if Grid.seq[i] != Grid.seq[i-1]
      return get_Q!(Grid, i) - get_Q!(Grid, i - 1)
    else
      return Δ{q}(Dict{Int64,Float64}(), Dict{Int64,Float64}())
    end
  elseif i == 1
    return Δ{q}(get_Q!(Grid, i).n, get_Q!(Grid, i).w)
  end
end
function get_Δ!(Grid::NestedGrid{p,q}, i::Int) where {p, q <: NestedQuadratureRule}
  get!(() -> calc_Δ!(Grid, i), Grid.Δs, i)
end
function get_Δ_weight!(Grid::NestedGrid, i::Int, j::Int)
  get_Δ!(Grid, i).w[j]
end
function get_Δ_weight!(Grid::NestedGrid, sv::SVector{p,Int64}, tup::Tuple) where {p}
  out = 1.0
  for i ∈ 1:p
    out *= get_Δ!(Grid, sv[i]).w[tup[i]]
  end
  out
end
#get_Δprod takes a vector of rule-sizes
#Δprods themselves are dictionaries with keys that have ints for node location.
function get_Δ_prod!(Grid::NestedGrid{p,q}, i::SVector{p, Int64}) where {p,q<:NestedQuadratureRule}
  get!(() -> calc_Δ_prod!(Grid, i), Grid.Δ_prods, i)
end
function expand_grid!(Grid::NestedGrid{p,q}, Δ_prod::Δprod{p}) where {p,q}
  for i ∈ keys(Δ_prod.w)
    Grid.grid.w[i] = get(Grid.grid.w, i, 0.0) + Δ_prod.w[i]
  end
end
function expand_grid!(Grid::NestedGrid{p,q}, j::SVector{p, Int64}) where {p,q}
  for i ∈ keys(get_Δ_prod!(Grid, j).w)
    Grid.grid.w[i] = get(Grid.grid.w, i, 0.0) + Grid.Δ_prods[j].w[i]
  end
end

@generated function calc_Δ_prod!(Grid::NestedGrid{p,q}, arg_indices::SVector{p,Int64}) where {p,q<:NestedQuadratureRule}
  quote
    out = Δprod($p)
    @nloops $p i dim -> begin
      keys(get_Δ!(Grid, arg_indices[dim]).w)
    end begin
      (@nref $p out i) = get_Δ_weight!( Grid, arg_indices, (@ntuple $p i) )
    end
    expand_grid!(Grid, out)
    out
  end
end
@generated function calc_Δ_prod(Grid::NestedGrid{p,q}, arg_indices::SVector{p,Int64}) where {p,q<:NestedQuadratureRule}
  quote
    out = Δprod($p)
    @nloops $p i dim -> begin
      keys(get_Δ!(Grid, arg_indices[dim]).w)
    end begin
      (@nref $p out i) = get_Δ_weight!( Grid, arg_indices, (@ntuple $p i) )
    end
    out
  end
end

function Base.length(Grid::NestedGrid)
  length(Grid.grid)
end
