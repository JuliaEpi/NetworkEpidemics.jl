
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
        update_state!(state, a, k, em)
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

"""
    average(em, state_0; kwargs)

Compute the average trajectory of an `AbstractEpidemicModel` over multiple stochastic simulations from the initial state `state_0`.

Keyword arguments:
* `nsims=1000`: the number of simulations.
* `nbins:100`: the number of time steps at which the averages should be computed. Time steps are uniformly distributed between `0` and `tmax`.
* `tmax=100.0`: the maximum time of each simulation.
* `nmax=1000`: the maximum number of iterations of each simulation.
* `progressbar=progressbar(em)`: if set to `true`, a progress bar will be displayed while the algorithm is running.
"""
function average(em::AbstractEpidemicModel, state_0; nsims=1000, nbins=100, tmax=100.0, nmax=1000, progressbar=progressbar(em))
    ts = LinRange(0.0, tmax, nbins)
    sums = init_sums(em, nbins)
    if progressbar
        p = Progress(nsims, dt=1.0)
    end
    for n in 1:nsims
        t, states = gillespie(em, state_0, tmax=tmax, nmax=nmax, process_output=false)
        l = 1
        for k in 1:nbins
            while l < length(t) && t[l+1] < ts[k]
                l += 1
            end
            update_sums!(sums, k, l, states, em)
        end
        if progressbar
            ProgressMeter.next!(p; showvalues = [(:iteration,n)])
        end
    end
    return finalize_sums(sums, nsims, em)
end

progressbar(::AbstractEpidemicModel) = false

function init_sums(em::AbstractEpidemicModel, nbins)
    dims = output_elem_size(em)
    sums = [zeros(output_type(em), dims) for i in 1:nbins]
    return sums
end

function update_sums!(sums, k, l, states, em::AbstractEpidemicModel)
    sums[k] .+= states[l]
end

function finalize_sums(sums, nsims, em::AbstractEpidemicModel)
    return float(sums)./nsims
end
