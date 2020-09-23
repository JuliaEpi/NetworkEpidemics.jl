module NetworkEpidemics

using ProgressMeter
using LinearAlgebra
using OrdinaryDiffEq
using LightGraphs
using LightGraphs.SimpleGraphs
using StatsBase


include("utils\\categorical_tree.jl")
include("utils\\rand.jl")
include("utils\\unzip.jl")
include("utils\\laplacian.jl")
include("utils\\arrays.jl")

export normalized_laplacian

include("models\\AbstractEpidemicModel.jl")

export AbstractEpidemicModel

include("EpidemicEvent.jl")

export EpidemicEvent, EmptyEpidemicEvent
export CPSimpleInfectionEvent, CPSimpleRecoveryEvent, CPSimpleRemovalEvent

include("models\\models\\metapopulation.jl")
include("models\\models\\contact_process.jl")
include("models\\models\\metaplex.jl")


include("gillespie.jl")
include("meanfield.jl")

export gillespie, average, meanfield

include("models\\dynamics\\AbstractCompartimentalModel.jl")
include("models\\dynamics\\SI.jl")
include("models\\dynamics\\SIS.jl")
include("models\\dynamics\\SIR.jl")

export Metapopulation, ContactProcess, Metaplex
export AbstractCompartimentalModel
export SI, SIS, SIR

end # module
