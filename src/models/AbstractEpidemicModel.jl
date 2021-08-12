
"""
    AbstractEpidemicModel

An abstract type representing an epidemic model.
"""
abstract type AbstractEpidemicModel end


"""
    init_state(em, t)

Setup the initial state for the Gillespie simulation.
"""
function init_state(em::AbstractEpidemicModel, state_0)
    return copy(state_0)
end

"""
    rate(em::AbstractEpidemicModel, k, state)

Return the reaction rate for reaction `k` for the epidemic model `em` at state `state`.
"""
function rate end

"""
    init_rates(em, state)

Compute the rates vector `a` for the state `state`.
"""
function init_rates end

"""
    init_output(em, state, n)

Initialize the output of the algorithm, allocating it and processing the initial `state`.
"""
function init_output(em::AbstractEpidemicModel, state::T, nmax) where T
    output = Vector{T}(undef, nmax)
    output[1] = copy(state)
    return output
end

"""
    update_state!(em, state, n)

Update the `state` of the simulation and the rates vector `a` with reaction `k` from the rates vector.
"""
function update_state_and_rates! end

"""
    update_output!(output, state, n, k, em)

Update the `output` of the simulation.

Arguments:
* `output`: the output array. `output[n]` gets mutated.
* `state`: the current state of the simulation
* `n`: the current time step.
* `k`: ...
"""
function update_output!(output, state, n, k, em::AbstractEpidemicModel)
    # NB. This is inefficient, so you'll want to overload it in most cases
    output[n] = copy(state) 
end

"""
    finalize_output(em, output)

Finalize the output.

Overwrite this to massage your output for easier plotting.
"""
finalize_output(em::AbstractEpidemicModel, output) = output


