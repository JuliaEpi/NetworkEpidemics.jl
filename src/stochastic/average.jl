
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
            while l < length(t) && t[l+1] <= ts[k]
                l += 1
            end
            update_sums!(sums, k, l, states, em)
        end
        if progressbar
            ProgressMeter.next!(p; showvalues = [(:iteration,n)])
        end
    end
    return ts, finalize_sums(sums, nsims, em)
end

"""
    progressbar(::AbstractEpidemicModel)

 Return whether a progress bar should be displayed by default when using `average`.
"""
progressbar(::AbstractEpidemicModel) = false

"""
    init_sums(em, nbins)

Return an array of size `nbins` of `em`'s output type, initialised to 0.
"""
function init_sums(em::AbstractEpidemicModel, nbins)
    dims = output_elem_size(em)
    sums = [zeros(output_type(em), dims) for i in 1:nbins]
    return sums
end

"""
    update_sums!(sums, k, l, states, em)

Add the information from `states[l]` to `sums[k]`.
"""
function update_sums!(sums, k, l, states, em::AbstractEpidemicModel)
    sums[k] .+= states[l]
end

"""
    finalize_sums(sums, nsims, em)

Perform some final operations on the agregated data from `sums`, e.g. dividing by `nsims` to get the mean.
"""
function finalize_sums(sums, nsims, em::AbstractEpidemicModel)
    return float.(finalize_output(em, sums))./nsims
end
