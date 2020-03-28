module NetworkEpidemics

using ProgressMeter
using LinearAlgebra
using DifferentialEquations
using LightGraphs
using LightGraphs.SimpleGraphs


include("gillespie.jl")
include("meanfield.jl")
include("utils\\rand.jl")
include("utils\\categorical_tree.jl")

include("models\\AbstractEpidemicModel.jl")
include("models\\models\\metapopulation.jl")
include("models\\models\\metaplex.jl")
include("models\\models\\contact_process.jl")

include("models\\dynamics\\SI.jl")
include("models\\dynamics\\SIS.jl")
include("models\\dynamics\\SIR.jl")

export AbstractEpidemicModel
export gillespie, meanfield
export Metapopulation, ContactProcess, Metaplex
export SI, SIS, SIR

end # module
