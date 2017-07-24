struct Δprod{p}
  w::Dict{SVector{p, Int64},Float64}
end
function Δprod(p)
  Δprod{p}(Dict{SVector{p, Int64},Float64}())
end
function Base.:+{p}(x::Δprod{p}, y::Δprod{p})
  Δprod(x.w + y.w)
end
function Base.:+{p}(x::Δprod{p}...)
  Δprod(sum(x.w))
end
function Base.setindex!{p}(A::Δprod{p}, x::Float64, i::Int64...)
  A.w[SVector{p,Int64}(i)] = x
end
function Base.getindex{p}(A::Δprod{p}, i::SVector{p,Int64})
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
function NestedGrid(p::Int64, q::DataType)
  NestedGrid{p, q}(Dict{Int64,q}(), Dict{Int64,Δ{q}}(), Dict{SVector{p, Int64}, Δprod}(), Dict{Int64,Float64}(), default(q), Δprod(p))
end
function NestedGrid(p::Int64, q::DataType, seq::Array{Int64,1})
  NestedGrid{p, q}(Dict{Int64,q}(), Dict{Int64,Δ{q}}(), Dict{SVector{p, Int64}, Δprod}(), Dict{Int64,Float64}(), seq, Δprod(p))
end
function get_node(Grid::NestedGrid, i::Int)
  Grid.node_value[i]
end
function get_Q!(Grid::NestedGrid{p, q}, i::Int) where {p, q <: NestedQuadratureRule}
  try
    return Grid.Qs[i]
  catch
    Grid.Qs[i] = q(Grid.seq[i])
    if length(Grid.Qs[i].n) > length(Grid.node_value)
      merge!(Grid.node_value, Grid.Qs[i].n)
    end
    return Grid.Qs[i]
  end
end
function get_Δ!{p, q <: NestedQuadratureRule}(Grid::NestedGrid{p,q}, i::Int)#
  try
    return Grid.Δs[i]
  catch
    if i > 1
      if Grid.seq[i] != Grid.seq[i-1]
        return Grid.Δs[i] = get_Q!(Grid, i) - get_Q!(Grid, i - 1)
      else
        return Grid.Δs[i] = Δ{q}(Dict{Int64,Float64}(), Dict{Int64,Float64}())
      end
    elseif i == 1
      return Grid.Δs[i] = Δ{q}(get_Q!(Grid, i).n, get_Q!(Grid, i).w)
    end
  end
end
function get_Δ_weight!(Grid::NestedGrid, i::Int, j::Int)
  get_Δ!(Grid, i).w[j]
end
function get_Δ_weight!{p}(Grid::NestedGrid, sv::SVector{p,Int64}, tup::Tuple)
  out = 1.0
  for i ∈ 1:p
    out *= get_Δ!(Grid, sv[i]).w[tup[i]]
  end
  out
end
#get_Δprod takes a vector of rule-sizes
#Δprods themselves are dictionaries with keys that have ints for node location.
function get_Δ_prod!{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q}, i::SVector{p, Int64})
  try
    return Grid.Δ_prods[i]
  catch
    return Grid.Δ_prods[i] = calc_Δ_prod!(Grid, i)
  end
end


@generated function calc_Δ_prod!{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q}, arg_indices::SVector{p,Int64})
  quote
    out = Δprod($p)
    @nloops $p i dim -> begin
      keys(get_Δ!(Grid, arg_indices[dim]).w)
    end begin
      (@nref $p out i) = get_Δ_weight!( Grid, arg_indices, (@ntuple $p i) )
    end
    Grid.Δ_prods[arg_indices] = out
    for i ∈ keys(out.w)
      Grid.grid.w[i] = get(Grid.grid.w, i, 0.0) + out.w[i]
    end
  end
end

function Base.length(Grid::NestedGrid)
  length(Grid.grid)
end
