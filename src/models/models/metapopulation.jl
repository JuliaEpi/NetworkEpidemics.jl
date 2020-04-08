
struct Metapopulation{T} <: AbstractEpidemicModel where T
    h::AbstractSimpleGraph # Population network
    D::Vector{Float64}
    dynamics::T
end

# handy little wrapper
Metapopulation(h, D::Real, dynamics) = Metapopulation(
    h,
    fill(float(D),
    num_states(dynamics)),
    dynamics
    )

# Implement the AbstractEpidemicModel interface

function init_rates(mp::Metapopulation, state)
    N = nv(mp.h)
    a =  zeros(N) # Top level reaction rates per Population node
    init_rates!(a, state, mp) # pass it down to the dynamics
    return a
end

output_type(::Metapopulation) = Int


# Not part of the interface
