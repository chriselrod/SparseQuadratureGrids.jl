@generated function smolyak!{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q}, l::Int)
  quote
    eval(parse("j_"*string($p+1)*" = "*string(l)))
    eval(parse("s_"*string($p+1)*" = 0"))
    @nloops $p i d -> begin
      1:j_{d+1}
    end d -> begin
      s_d = s_{d+1} + i_d - 1
      j_d = l - s_d
    end begin
      calc_Δ_prod!(Grid, SVector((@ntuple $p i)))
    end
  end
end

@generated function l2!{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q}, l::Int)
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
      calc_Δ_prod!(Grid, SVector((@ntuple $p i)))
    end
  end
end

@generated function lx!{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q}, l::Int, x::Real)
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
      calc_Δ_prod!(Grid, SVector((@ntuple $p i)))
    end
  end
end
