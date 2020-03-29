
struct Metapopulation{T} <: AbstractEpidemicModel where T
    h::AbstractSimpleGraph # Population network
    D::Union{Float64, Vector{Float64}}
    dynamics::T
end

# Implement the AbstractEpidemicModel interface

function init_state(mp::Metapopulation, state_0, t)
    state = copy(state_0)
    return state
end

function init_rates(mp::Metapopulation, state)
    N = nv(mp.g)
    a = Vector{Float64}(undef, N) # Top level reaction rates per Population node
    init_rates!(a, state, mp) # pass it down to the dynamics
    return a
end

#update_state!(state, a, mp::Metapopulation, k) = _NI("update_state!")

# Not part of the interface
