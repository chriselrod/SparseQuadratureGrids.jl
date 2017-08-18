module SparseQuadratureGrids

using JLD2
using Base.Cartesian
using StaticArrays

import Base:    length,
                getindex,
                setindex!,
                Val,
                convert

export  SplitWeights,
        QuadratureRule,
        KronrodPatterson,
        GenzKeister,
        NestedGrid,
        calc_Δ_prod!,
        smolyak!,
        FlatGrid,
        FlattenedGrid,
        GridContainer,
        GridBuild,
        AdaptiveBuild,
        aPrioriBuild,
        RawBuild,
        CacheBuild,
        Adaptive,
        AdaptiveRaw,
        Smolyak,
        SmolyakRaw,
        eval_grid!,
        calc_grid!,
        GridVessel,
        DynamicGridVessel,
        StaticGridVessel,
        default,
        SlidingVecFun,
        SlidingVec,
        SlidingVector,
        MatrixVecSVec,
        set!,
        eval!

#BLAS.set_num_threads(1)
include("univariate_rules.jl")
include("sliding_vectors.jl")
include("grid.jl")
include("grid_builders.jl")
include("grid_containers.jl")
include("grid_interfacing.jl")

end # module
