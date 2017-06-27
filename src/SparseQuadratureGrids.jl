module SparseQuadratureGrids

using JLD
using Base.Cartesian
using StaticArrays

import Base.length

export  SplitWeights,
        QuadratureRule,
        KronrodPatterson,
        GenzKeister,
        NestedGrid,
        calc_Δ_prod!,
        smolyak!

include("univariate_rules.jl")
include("grid.jl")
include("grid_builders.jl")
include("grid_interfacing.jl")

end # module
