module NetworkEpidemics

using ProgressMeter
using LinearAlgebra
using OrdinaryDiffEq
using LightGraphs
using LightGraphs.SimpleGraphs
using StatsBase
using DataStructures
using FillArrays


include("utils/categorical_tree.jl")
include("utils/rand.jl")
include("utils/unzip.jl")
include("utils/laplacian.jl")
include("utils/arrays.jl")
include("utils/NotImplementedError.jl")

export normalized_laplacian

include("models/AbstractEpidemicModel.jl")

export AbstractEpidemicModel

include("models/models/metapopulation.jl")
include("models/models/contact_process.jl")
include("models/models/metaplex.jl")
include("models/models/heterogeneous_metaplex.jl")


include("stochastic/AbstractSimulator.jl")
include("stochastic/gillespie.jl")
include("stochastic/average.jl")
include("deterministic/meanfield.jl")

export gillespie, average, meanfield

include("models/dynamics/AbstractCompartimentalModel.jl")
include("models/dynamics/SI.jl")
include("models/dynamics/SIS.jl")
include("models/dynamics/SIR.jl")

export Metapopulation, ContactProcess, Metaplex
export HeterogeneousMetaplex, HeterogeneousMetapopulation
export AbstractCompartimentalModel
export SI, SIS, SIR

end # module
