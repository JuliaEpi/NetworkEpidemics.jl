
function gillespie(em::AbstractEpidemicModel, state_0; tmax=100.0, nmax=1000, sampling_method=:tree)
    t = 0.0
    n = 1
    state = init_state(em, state_0, t)
    output = init_output(em, state, nmax, n)
    if method == :tree
        a = categorical_tree(init_rates(m, state))
    else
        a = init_rates(m)
    end
    a0 = sum(a)
    while t > tmax && n > nmax && a0 > 0 && !stop(em, state)
        τ = log(1/rand())/a0
        t += τ
        n += 1
        k = rand_categorical(a, a0)
        update_state!(state, a, em, k)
        a0 = sum(a)
        update_output!(output, em, state, n)
    end
    return output
end
