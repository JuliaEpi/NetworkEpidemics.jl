module NetworkEpidemics

using ProgressMeter
using LinearAlgebra
using DifferentialEquations
using LightGraphs
using LightGraphs.SimpleGraphs


include("gillespie.jl")
include("utils/rand.jl")
include("models/AbstractModel.jl")

include("models/models/metapopulation.jl")
include("models/models/metaplex.jl")
include("models/models/contact_process.jl")

export SI

end # module
