# NetworkEpidemics.jl
A Julia package for simulating epidemics on various network models based on [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl).

This currently includes simulation of SI, SIS and SIR dynamics on Metapopulation, Contact Processes and a brand new model named Metaplex (name may be subject to change).

Each model can be simulated using either a generic Gillespie algorithm or by integrating meanfield equations.

## Installation
```julia
pkg> add https://github.com/csimal/NetworkEpidemics.jl
```

## Planned features
* An easy way to specify an arbitrary compartimental model that can then directly be simulated.
* Support for [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) (Related to previous point).
* Support for more deterministic models, e.g. the Degree-Based approximation.
* Parallel computation of stochastic averages.