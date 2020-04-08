
"""
    Metaplex{T} <: AbstractEpidemicModel

A type representing a metaplex, i.e. a metapopulation model with an embedded individual network, with arbitrary dynamics.
"""
struct Metaplex{T} <: AbstractEpidemicModel where T
    g::AbstractSimpleGraph # individual relation network
    h::AbstractSimpleGraph # population network
    D::Vector{Float64}
    dynamics::T
end

Metaplex(g, h, D::Real, dynamics) = Metaplex(
    g,
    h,
    fill(float(D), num_states(dynamics)),
    dynamics
    )

Metapopulation(mpx::Metaplex) = Metapopulation(mpx.h, mpx.D, mpx.dynamics)

ContactProcess(mpx::Metaplex) = ContactProcess(mpx.g, mpx.dynamics)

function init_rates(mpx::Metaplex, state)
    N = nv(mpx.g)
    a = zeros(N)
    init_rates!(a, state, mpx)
    return a
end

output_type(::Metaplex) = Int
