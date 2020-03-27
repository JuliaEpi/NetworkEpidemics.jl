
struct Metapopulation{T} <: AbstractEpidemicModel
    h::AbstractSimpleGraph # Population network
    dynamics::T
end

function init_state(mp::Metapopulation, state_0, t)

end

function init_rates(mp::Metapopulation, state)

end

function init_output(mp::Metapopulation, state, n)

end

function update_state!(state, a, mp::Metapopulation, k)

end

function update_output!(output, state, mp::Metapopulation, n)

end
