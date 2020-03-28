
struct SI
    β::Real
end

rate(i, β, D, state, mp) = β*state[i,1]*state[i,2] + sum( D .* state[i,:] ) * outdegree(mp.g, i)

function init_rates!(a, state, mp::Metapopulation{SI})
    N = nv(mp.g)
    D = mp.D
    β = mp.dynamics.β
    for i in 1:N
        a[i] = rate(i, β, D, state, mp)
    end
end

function update_state!(state, a, k, mp::Metapopulation{SI})
    N = nv(mp.g)
    D = mp.D
    β = mp.dynamics.β
    p = vcat(β*state[k,1]*state[k,2], D.* state[k,:] * outdegree(mp.g, k))
    j = sample(ProbabilityWeights(p))
    if j == 1
        state[k,1] -= 1
        state[k,2] += 1
        a[k] = rate(k, β, D, state, mp)
    else # migration from node μ to node ν
        μ = k
        ν = rand(outneighbors(mp.g, μ))
        state[μ,j-1] -= 1
        state[ν,j-1] += 1
        a[μ] = rate(μ, β, D, state, mp)
        a[ν] = rate(ν, β, D, state, mp)
    end
end
