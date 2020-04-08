
"""
    SI

A type representing a Susceptible Infected model with infection parameter `β`.
"""
struct SI <: AbstractCompartimentalModel
    β::Real
end

num_states(::SI) = 2

# Metapopulation

function rate(i, state, mp::Metapopulation{SI})
    β = mp.dynamics.β
    D = mp.D
    β*state[i,1]*state[i,2] + sum( D .* state[i,:] ) * outdegree(mp.h, i)
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
        a[k] = rate(k, state, mp)
    else # migration from node μ to node ν
        μ = k
        ν = rand(outneighbors(mp.h, μ))
        state[μ,j-1] -= 1
        state[ν,j-1] += 1
        a[μ] = rate(μ, state, mp)
        a[ν] = rate(ν, state, mp)
    end
end

function meanfield_fun(mp::Metapopulation{SI})
    L = normalized_laplacian(mp.h)
    β = mp.dynamics.β
    D = mp.D
    f! = function(dx, x, p, t)
        tmp = β.*x[:,1].*x[:,2]
        dx[:,1] .= -tmp - D[1]*L'*x[:,1]
        dx[:,2] .= tmp - D[2]*L'x[:,2]
    end
    return f!
end

# Contact Process

function rate(k, state, cp::ContactProcess{SI})
    if !state[k]
        l = length(filter(i->state[i], inneighbors(cp.g, k)))
        return cp.dynamics.β * l
    else
        return 0.0
    end
end

function init_output(cp::ContactProcess{SI}, state, nmax)
    output = Vector{Tuple{Int,EpidemicEvent}}(undef, nmax)
    output[1] = (sum(state), EmptyEpidemicEvent())
    return output
end

function update_state!(state, a, k, cp::ContactProcess{SI})
    β = cp.dynamics.β
    state[k] = true
    a[k] = 0.0
    for i in Iterators.filter(j->!state[j], outneighbors(cp.g, k))
        a[i] += β
    end
end

function update_output!(output, state, n, k, cp::ContactProcess{SI})
    old = output[n-1][1]
    # pick a random infected neighbor as the one who infected k
    j = rand(filter(i->state[i], inneighbors(cp.g, k)))
    output[n] = (old+1, CPSimpleInfectionEvent(j,k))
end

function meanfield_fun(cp::ContactProcess{SI})
    A = adjacency_matrix(cp.g)
    β = cp.dynamics.β
    f! = function(dx, x, p, t)
        dx .= β.*(1.0 .- x).(A*x)
    end
    return f!
end

# Metaplex

function rate(k, state, mpx::Metaplex{SI})
    β : mpx.dynamics.β
    D = mpx.D
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    popcounts = state[:popcounts]
    if !x_i[k]
        l = length(filter(i->x_i[i] && x_μ[i]==x_μ[k], inneighbors(mpx.g, k)))
        return β*l + D[1]
    else
        return D[2]
    end
end

function update_state!(state, a, k, mpx::Metaplex{SI})
    N = nv(mpx.g)
    M = nv(mpx.h)
    g = mpx.g
    h = mpx.h
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    popcounts = state[:popcounts]
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

function meanfield_fun(mpx::Metaplex{SI})
    N = nv(mpx.g)
    M = nv(mpx.h)
    A = adjacency_matrix(mpx.g)
    L = normalized_laplacian(mpx.h)
    β = mpx.dynamics.β
    D = mpx.D
    f! = function(dx, x, p, t)
        for i in 1:N
            dx[1,i,:] .= .-D[1]*L'*x[1,i,:]
            dx[2,i,:] .= .-D[2]*L'*x[2,i,:]
        end
        for μ in 1:M
            infection = β*(x[1,:,μ].*(A*u[2,:,μ]))
            dx[1,:,μ] .-= infection
            dx[2,:,μ] .== infection
        end
    end
    return f!
end
