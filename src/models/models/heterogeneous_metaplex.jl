
"""
    HeterogeneousMetaplex{T} <: AbstractEpidemicModel

A type representing a metaplex, i.e. a metapopulation model with embedded individual networks, with arbitrary dynamics.
"""
struct HeterogeneousMetaplex{T} <: AbstractEpidemicModel
    g::Vector{<:AbstractGraph} # individual relation network
    h::AbstractSimpleGraph # population network
    D::Vector{Float64}
    dynamics::T
end

HeterogeneousMetaplex(g, h, D::Real, dynamics) = HeterogeneousMetaplex(
    g,
    h,
    fill(float(D), num_states(dynamics)),
    dynamics
    )

Metapopulation(mpx::HeterogeneousMetaplex) = Metapopulation(mpx.h, mpx.D, mpx.dynamics)

#ContactProcess(mpx::HeterogeneousMetaplex) = ContactProcess(mpx.g, mpx.dynamics)

function init_rates(mpx::HeterogeneousMetaplex, state)
    N = nv(mpx.g[1])
    a = zeros(N)
    init_rates!(a, state, mpx)
    return a
end

output_type(::HeterogeneousMetaplex) = Int

progressbar(::HeterogeneousMetaplex) = true


struct HeterogeneousMetapopulation{T} <: AbstractEpidemicModel
    h::AbstractSimpleGraph
    N::Int
    ks::Vector{Int}
    D::Vector{Float64}
    dynamics::T
end

HeterogeneousMetapopulation(h, N, ks, D::Real, dynamics) = HeterogeneousMetapopulation(h, N, ks, fill(float(D), num_states(dynamics)), dynamics)

HeterogeneousMetapopulation(h, N, ks::Int, D, dynamics) = HeterogeneousMetapopulation(h, N, fill(ks, nv(h)), D, dynamics)

function init_rates(mp::HeterogeneousMetapopulation, state)
    a = zeros(nv(mp.h))
    init_rates!(a, state, mp)
    return a
end

output_type(::HeterogeneousMetapopulation) = Int

