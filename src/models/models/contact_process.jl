
"""
    ContactProcess{T} <: AbstractEpidemicModel

A type representing a contact process on network `g` with arbitrary dynamics.
"""
struct ContactProcess{T} <: AbstractEpidemicModel where T
    g::AbstractSimpleGraph
    dynamics::T
end

function init_rates(cp::ContactProcess, state)
    N = nv(cp.g)
    a = Vector{Float64}(undef, N)
    init_rates!(a, state, cp)
    return a
end
