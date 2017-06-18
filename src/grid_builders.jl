@generated function smolyak!{p,q<:NestedQuadratureRule}(Grid::NestedGrid{p,q}, l::Int64)
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
