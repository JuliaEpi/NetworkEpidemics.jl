
"""
    SI

A type representing a Susceptible Infected model with infection parameter `β`.
"""
struct SI
    β::Real
end

num_states(::SI) = 2

# Metapopulation

_rate(i, β, D, state, mp::Metapopulation{SI}) = β*state[i,1]*state[i,2] + sum( D .* state[i,:] ) * outdegree(mp.h, i)

function init_rates!(a, state, mp::Metapopulation{SI})
    N = nv(mp.h)
    D = mp.D
    β = mp.dynamics.β
    for i in 1:N
        a[i] = _rate(i, β, D, state, mp)
    end
end

function init_output(mp::Metapopulation{SI}, state, nmax)
    N = nv(mp.h)
    output = [Array{Int,2}(undef, N, 2) for i in 1:nmax]
    output[1] .= state
    return output
end

function update_state!(state, a, k, mp::Metapopulation{SI})
    N = nv(mp.h)
    D = mp.D
    β = mp.dynamics.β
    p = vcat(β*state[k,1]*state[k,2], D.* state[k,:] * outdegree(mp.h, k))
    j = sample(ProbabilityWeights(p))
    if j == 1
        state[k,1] -= 1
        state[k,2] += 1
        a[k] = _rate(k, β, D, state, mp)
    else # migration from node μ to node ν
        μ = k
        ν = rand(outneighbors(mp.h, μ))
        state[μ,j-1] -= 1
        state[ν,j-1] += 1
        a[μ] = _rate(μ, β, D, state, mp)
        a[ν] = _rate(ν, β, D, state, mp)
    end
end

function update_output!(output, state, n, k, mp::Metapopulation{SI})
    output[n] .= state
end

function finalize_output(mp::Metapopulation{SI}, output)
    N = nv(mp.h)
    n = length(output)
    s = Array{Int,2}(undef, n, N)
    i = Array{Int,2}(undef, n, N)
    for k in 1:n
        s[k,:] .= output[k][:,1]
        i[k,:] .= output[k][:,2]
    end
    return s, i
end

# Contact Process

function init_rates!(a, state, cp::ContactProcess{SI})
    N = nv(cp.g)
    for k in 1:N
        if !state[k]
            l = length(filter(i->state[i], inneighbors(cp.g, k)))
            a[k] = cp.dynamics.β * l
        end
    end
end

function init_output(cp::ContactProcess{SI}, state, nmax)
    output = Vector{Tuple{Int,EpidemicEvent}}(undef, nmax)
    output[1] = (sum(state), EmptyEpidemicEvent())
    return output
end

function update_state!(state, a, k, cp::ContactProcess{SI})
    state[k] = true
    a[k] = 0.0
    for i in Iterators.filter(j->!state[j], outneighbors(g, k))
        a[i] += β
    end
end

function update_output!(output, state, n, k, cp::ContactProcess{SI})
    # pick a random infected neighbor as the one who infected k
    j = rand(filter(i->state[i], inneighbors(cp.g, k)))
    output[n] = (n, CPSimpleInfectionEvent(j,k))
end

function finalize_output(cp::ContactProcess{SI}, output)
    return unzip(output)
end

# Metaplex

function init_state(mpx::Metaplex{SI}, state_0, t)
    N = nv(mpx.g)
    M = nv(mpx.h)
    state_i = copy(state_0[1]) # epidemic state
    state_μ = copy(state_0[2]) # location
    popcounts = zeros(Int, 2, M)
    for k in 1:N
        if !state_i[k]
            popcounts[1, state_μ[k]] += 1
        else
            popcounts[2, state_μ[k]] += 1
        end
    end
    return state_i, state_μ, popcounts
end

function init_rates!(a, state, mpx::Metaplex{SI})
    N = nv(mpx.g)
    x_i, x_μ, popcounts = state
    β = mpx.dynamics.β
    D = mpx.D
    for k in 1:N
        if !state_i[k]
            l = length(filter(i->x_i[i] && x_μ[i] == x_μ[k], inneighbors(mpx.g, k)))
            a[k] = β*l + D[1]
        else
            a[k] = D[2]
        end
    end
end

function init_output(mpx::Metaplex{SI}, state, nmax)
    M = nv(mpx.h)
    output = [Array{Int,2}(undef, 2, M) for i in 1:nmax]
    output[1] .= state[3]
    return output
end

function update_state!(state, a, k, mpx::Metaplex{SI})
    N = nv(mpx.g)
    M = nv(mpx.h)
    g = mpx.g
    h = mpx.h
    x_i, x_μ, popcounts = state
    β = mpx.dynamics.β
    D = mpx.D
    μ = x_μ[k]
    if !x_i[k]
        if rand()*a[k] < D[1]
            ν = rand(outneighbors(h,μ))
            x_μ[k] = ν
            popcounts[1,μ] -= 1
            popcounts[1,ν] += 1
            l = length(filter(i->x_i[i] && x_μ[i]==ν, inneighbors(g,k)))
            a[k] = β*l + D[1]
        else
            x_i[k] = true
            popcounts[1,μ] -= 1
            popcounts[2,μ] += 1
            a[k] = D[2]
            for i in Iterators.filter(j->!x_i[j] && x_μ[j]==μ, outneighbors(g,k))
                a[i] += β
            end
        end
    else
        ν = rand(outneighbors(h,μ))
        x_μ[k] = ν
        popcounts[2,μ] -= 1
        popcounts[2,ν] += 1
        for i in Iterators.filter(j->!x_i[j], outneighbors(g,k))
            if x_μ[i] == μ
                a[i] -= β
            elseif x_μ[i] == ν
                a[i] += β
            end
        end
    end
end

function update_output!(output, state, n, k, mpx::Metaplex{SI})
    output[n] .= state[3]
end
