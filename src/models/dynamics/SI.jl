
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
    β*state[i,1]*state[i,2] + sum( D .* state[i,:] ) * (outdegree(mp.h, i) > 0)
end

function update_state_and_rates!(state, a, k, mp::Metapopulation{SI})
    N = nv(mp.h)
    D = mp.D
    β = mp.dynamics.β
    p = vcat(β*state[k,1]*state[k,2], D.* state[k,:] * (outdegree(mp.h, k) > 0))
    j = sample(ProbabilityWeights(p))
    if j == 1 # infection in node k
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
        dx[:,2] .= tmp - D[2]*L'*x[:,2]
    end
    return f!
end

# Contact Process

function rate(k, state, cp::ContactProcess{SI})
    if state[k] == 1
        l = length(filter(i->state[i]==2, inneighbors(cp.g, k)))
        return cp.dynamics.β * l
    else
        return 0.0
    end
end

function update_state_and_rates!(state, a, k, cp::ContactProcess{SI})
    β = cp.dynamics.β
    state[k] = 2
    a[k] = 0.0
    for i in Iterators.filter(j->state[j]==1, outneighbors(cp.g, k))
        a[i] += β
    end
end

function update_output!(output, state, n, k, cp::ContactProcess{SI})
    output[n] = output[n-1] + [-1, 1]
end

function init_state_mf(cp::ContactProcess{SI}, x0)
    return float(x0 .- 1)
end

function meanfield_fun(cp::ContactProcess{SI})
    A = adjacency_matrix(cp.g)
    β = cp.dynamics.β
    f! = function(dx, x, p, t)
        dx .= β.*(1.0 .- x).*(A*x)
    end
    return f!
end

function finalize_meanfield(cp::ContactProcess{SI}, sol)
    N = nv(cp.g)
    n = length(sol.u)
    output = Array{Float64,2}(undef, n, N)
    for i in 1:n
        output[i,:] .= sol.u[i]
    end
    return sol.t, output
end

# Metaplex

function rate(k, state, mpx::Metaplex{SI})
    β = mpx.dynamics.β
    D = mpx.D
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    #popcounts = state[:popcounts]
    if x_i[k] == 1
        l = length(filter(i->x_i[i]==2 && x_μ[i]==x_μ[k], inneighbors(mpx.g, k)))
        return β*l + D[1]
    else
        return D[2]
    end
end

function update_state_and_rates!(state, a, k, mpx::Metaplex{SI})
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
    if x_i[k] == 1
        if rand()*a[k] < D[1]
            ν = rand(outneighbors(h,μ))
            x_μ[k] = ν
            popcounts[μ,1] -= 1
            popcounts[ν,1] += 1
            l = length(filter(i->x_i[i]==2 && x_μ[i]==ν, inneighbors(g,k)))
            a[k] = β*l + D[1]
        else
            x_i[k] = 2
            popcounts[μ,1] -= 1
            popcounts[μ,2] += 1
            a[k] = D[2]
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g,k))
                a[i] += β
            end
        end
    else
        ν = rand(outneighbors(h,μ))
        x_μ[k] = ν
        popcounts[μ,2] -= 1
        popcounts[ν,2] += 1
        for i in Iterators.filter(j->x_i[j]==1, outneighbors(g,k))
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
            infection = β*(x[1,:,μ].*(A*x[2,:,μ]))
            dx[1,:,μ] .-= infection
            dx[2,:,μ] .+= infection
        end
    end
    return f!
end

# Heterogeneous Metaplex

function rate(k, state, mpx::HeterogeneousMetaplex{SI})
    β = mpx.dynamics.β
    D = mpx.D
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    μ = x_μ[k]
    #popcounts = state[:popcounts]
    if x_i[k] == 1
        l = length(filter(i->x_i[i]==2 && x_μ[i]==μ, 
                inneighbors(mpx.g[μ], k))
            )
        return β*l + D[1]
    else
        return D[2]
    end
end

function update_state_and_rates!(state, a, k, mpx::HeterogeneousMetaplex{SI})
    g = mpx.g
    h = mpx.h
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    popcounts = state[:popcounts]
    β = mpx.dynamics.β
    D = mpx.D
    μ = x_μ[k]
    if x_i[k] == 1
        if rand()*a[k] < D[1]
            ν = rand(outneighbors(h,μ))
            x_μ[k] = ν
            popcounts[μ,1] -= 1
            popcounts[ν,1] += 1
            l = length(filter(i->x_i[i]==2 && x_μ[i]==ν, inneighbors(g[μ],k)))
            a[k] = β*l + D[1]
        else
            x_i[k] = 2
            popcounts[μ,1] -= 1
            popcounts[μ,2] += 1
            a[k] = D[2]
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g[μ],k))
                a[i] += β
            end
        end
    else
        ν = rand(outneighbors(h,μ))
        x_μ[k] = ν
        popcounts[μ,2] -= 1
        popcounts[ν,2] += 1
        for i in Iterators.filter(j->x_i[j]==1, outneighbors(g[μ],k))
            if x_μ[i] == μ
                a[i] -= β
            elseif x_μ[i] == ν
                a[i] += β
            end
        end
    end
end

function meanfield_fun(mpx::HeterogeneousMetaplex{SI})
    N = nv(mpx.g[1])
    M = nv(mpx.h)
    A = adjacency_matrix.(mpx.g)
    L = normalized_laplacian(mpx.h)
    β = mpx.dynamics.β
    D = mpx.D
    f! = function(dx, x, p, t)
        for i in 1:N
            dx[1,i,:] .= .-D[1]*L'*x[1,i,:]
            dx[2,i,:] .= .-D[2]*L'*x[2,i,:]
        end
        for μ in 1:M
            infection = β*(x[1,:,μ].*(A[μ]*x[2,:,μ]))
            dx[1,:,μ] .-= infection
            dx[2,:,μ] .+= infection
        end
    end
    return f!
end
