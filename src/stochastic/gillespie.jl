
"""
    gillespie(em::AbstractEpidemicModel, state_0; kwargs)

Compute an exact trajectory of an `AbstractEpidemicModel` with initial state `state_0` using a generic gillespie algorithm.

`em` must implement the `AbstractEpidemicModel` interface.
Keyword arguments:
* `tmax=100.0`: maximum time of the simulation
* `nmax=1000`: maximum number of iterations of the simulation
* `sampling_method=:tree`: Sampling method for choosing what reaction happens at each iteration. The default is `:tree`, which takes ``O(\\log N)`` per iteration, where ``N`` is the length of the reaction rate vector. Set to `:array` for a method with ``O(N)`` sampling time, but faster for small numbers of reactions (≈ under 1000).
"""
function gillespie(em::AbstractEpidemicModel, state_0; tmax=100.0, nmax=1000, sampling_method=:tree, should_stop=default_stop, process_output=true)
    t = 0.0
    n = 1
    state = init_state(em, state_0)
    output = init_output(em, state, nmax)
    ts = Vector{Float64}(undef, nmax)
    ts[1] = t
    if sampling_method == :tree
        a = CategoricalTree(init_rates(em, state))
    else
        a = init_rates(m, state)
    end
    a0 = sum(a)
    while t < tmax && n < nmax && a0 > 0 && !should_stop(em, state)
        τ = log(1/rand())/a0
        t += τ
        n += 1
        k = rand_categorical(a, a0)
        update_state_and_rates!(state, a, k, em)
        a0 = sum(a)
        update_output!(output, state, n, k, em)
        ts[n] = t
    end
    if process_output
        return ts[1:n], finalize_output(em, output[1:n])
    else
        return ts[1:n], output[1:n]
    end
end