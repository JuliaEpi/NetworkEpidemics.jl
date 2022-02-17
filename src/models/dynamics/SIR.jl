
"""
    SIR

A type representing a Susceptible-Infected-Recovered model with infection rate `β` and recovery rate `δ`.
"""
struct SIR <: AbstractCompartimentalModel
    β::Real
    δ::Real
end

SI(sir::SIR) = SI(sir.β)

num_states(::SIR) = 3

# Metapopulation

function rate(i, state, mp::Metapopulation{SIR})
    β = mp.dynamics.β
    δ = mp.dynamics.δ
    D = mp.D
    h = mp.h
    β*state[i,1]*state[i,2] + δ*state[i,2] + sum(D .* state[i,:]) * (outdegree(h,i)>0)
end

function update_state_and_rates!(state, a, k, mp::Metapopulation{SIR})
    M = nv(mp.h)
    D = mp.D
    β = mp.dynamics.β
    δ = mp.dynamics.δ
    h = mp.h
    p = vcat(β*state[k,1]*state[k,2], δ*state[k,2], D.*state[k,:]*(outdegree(h,k)))
    j = sample(ProbabilityWeights(p, a[k]))
    if j == 1 # infection
        state[k,1] -= 1
        state[k,2] += 1
        a[k] = rate(k, state, mp)
    elseif j == 2 # recovery
        state[k,2] -= 1
        state[k,3] += 1
        a[k] = rate(k, state, mp)
    else # migration
        i = rand(outneighbors(h, k))
        state[k,j-2] -= 1
        state[i,j-2] += 1
        a[k] = rate(k, state, mp)
        a[i] = rate(i, state, mp)
    end
end

function meanfield_fun(mp::Metapopulation{SIR})
    L = laplacian_matrix(mp.h)
    β = mp.dynamics.β
    δ = mp.dynamics.δ
    D = mp.D
    f! = function(dx, x, p, t)
        infection = β*x[:,1].*x[:,2]
        recovery = δ*x[:,2]
        dx[:,1] .= .-infection .- D[1].*(L'*x[:,1])
        dx[:,2] .= infection .- recovery .- D[2].*(L'*x[:,2])
        dx[:,3] .= recovery .- D[3].*(L'x[:,3])
    end
    return f!
end

# Contact Process

function rate(k, state, cp::ContactProcess{SIR})
    if state[k] == 1
        l = length(filter(i->state[i]==2, inneighbors(cp.g, k)))
        return cp.dynamics.β * l
    elseif state[k] == 2
        return cp.dynamics.δ
    else
        return 0.0
    end
end

function update_state_and_rates!(state, a, k, cp::ContactProcess{SIR})
    g = cp.g
    β = cp.dynamics.β
    if state[k] == 1
        state[k] = 2
        for i in Iterators.filter(j->state[j]==1, outneighbors(g, k))
            a[i] += β
        end
    elseif state[k] == 2
        state[k] = 3
        for i in Iterators.filter(j->state[j]==1, outneighbors(g, k))
            a[i] -= β
        end
    end
    a[k] = rate(k, state, cp)
end

function update_output!(output, state, n, k, cp::ContactProcess{SIR})
    if state[k] == 2
        output[n] = output[n-1] + [-1, 1, 0]
    else
        output[n] = output[n-1] + [0, -1, 1]
    end
end

function meanfield_fun(cp::ContactProcess{SIR})
    A = adjacency_matrix(cp.g)
    β = cp.dynamics.β
    δ = cp.dynamics.δ
    f! = function(dx, x, p, t)
        infection = β.*x[:,1].*(A*x[:,2])
        dx[:,1] .= .- infection
        dx[:,2] .= infection .- δ.*x[:,2]
        dx[:,3] .= δ.*x[:,2]
    end
    return f!
end

# Metaplex

function rate(k, state, mpx::Metaplex{SIR})
    β = mpx.dynamics.β
    δ = mpx.dynamics.δ
    D = mpx.D
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    od = outdegree(mpx.h, x_μ[k])
    if x_i[k] == 1
        l = length(filter(i->x_i[i]==2 && x_μ[i]==x_μ[k], inneighbors(mpx.g, k)))
        return β*l + D[1]*od
    elseif x_i[k] == 2
        return δ + D[2]*od
    else
        return D[3]*od
    end
end

function update_state_and_rates!(state, a, k, mpx::Metaplex{SIR})
    N = nv(mpx.g)
    M = nv(mpx.h)
    g = mpx.g
    h = mpx.h
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    od = outdegree(mpx.h, x_μ[k])
    popcounts = state[:popcounts]
    β = mpx.dynamics.β
    δ = mpx.dynamics.δ
    D = mpx.D
    μ = x_μ[k]
    if x_i[k] == 1
        if rand()*a[k] < D[1]*od # susceptible node migrates
            ν = rand(outneighbors(h, μ))
            x_μ[k] = ν
            popcounts[μ,1] -= 1
            popcounts[ν,1] += 1
            a[k] = rate(k, state, mpx)
        else # infection
            x_i[k] = 2
            popcounts[μ,1] -= 1
            popcounts[μ,2] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g,k))
                a[i] += β
            end
        end
    elseif x_i[k] == 2
        if rand()*a[k] < D[2]*od # infected node migrates
            ν = rand(outneighbors(h,μ))
            x_μ[k] = ν
            popcounts[μ,2] -= 1
            popcounts[ν,2] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1, outneighbors(g,k))
                if x_μ[i] == μ
                    a[i] -= β
                elseif x_μ[i] == ν
                    a[i] += β
                end
            end
        else # recovery
            x_i[k] = 3
            popcounts[μ,2] -= 1
            popcounts[μ,3] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g,k))
                a[i] -= β
            end
        end
    else # recovered node migrates
        ν = rand(outneighbors(h,μ))
        x_μ[k] = ν
        popcounts[μ,3] -= 1
        popcounts[ν,3] += 1
        a[k] = rate(k, state, mpx)
    end
end

function meanfield_fun(mpx::Metaplex{SIR})
    N = nv(mpx.g)
    M = nv(mpx.h)
    A = adjacency_matrix(mpx.g)
    L = laplacian_matrix(mpx.h)
    β = mpx.dynamics.β
    δ = mpx.dynamics.δ
    D = mpx.D
    f! = function(dx, x, p, t)
        for i in 1:N
            dx[1,i,:] .= .-D[1]*L'*x[1,i,:]
            dx[2,i,:] .= .-D[2]*L'*x[2,i,:]
            dx[3,i,:] .= .-D[3]*L'*x[3,i,:]
        end
        for μ in 1:M
            infection = β.*x[1,:,μ].*(A*x[2,:,μ])
            recovery = δ.*x[2,:,μ]
            dx[1,:,μ] .-= infection
            dx[2,:,μ] .+= infection .- recovery
            dx[3,:,μ] .+= recovery
        end
    end
    return f!
end

# Heterogeneous Metaplex

function rate(k, state, mpx::HeterogeneousMetaplex{SIR})
    β = mpx.dynamics.β
    δ = mpx.dynamics.δ
    D = mpx.D
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    μ = x_μ[k]
    od = outdegree(mpx.h, x_μ[k])
    if x_i[k] == 1
        l = length(filter(i->x_i[i]==2 && x_μ[i]==μ, inneighbors(mpx.g[μ], k)))
        return β*l + D[1]*od
    elseif x_i[k] == 2
        return δ + D[2]*od
    else
        return D[3]*od
    end
end

function update_state_and_rates!(state, a, k, mpx::HeterogeneousMetaplex{SIR})
    g = mpx.g
    h = mpx.h
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    od = outdegree(mpx.h, x_μ[k])
    popcounts = state[:popcounts]
    β = mpx.dynamics.β
    D = mpx.D
    μ = x_μ[k]
    if x_i[k] == 1
        if rand()*a[k] < D[1]*od # susceptible node migrates
            ν = rand(outneighbors(h, μ))
            x_μ[k] = ν
            popcounts[μ,1] -= 1
            popcounts[ν,1] += 1
            a[k] = rate(k, state, mpx)
        else # infection
            x_i[k] = 2
            popcounts[μ,1] -= 1
            popcounts[μ,2] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g[μ],k))
                a[i] += β
            end
        end
    elseif x_i[k] == 2
        if rand()*a[k] < D[2]*od # infected node migrates
            ν = rand(outneighbors(h,μ))
            x_μ[k] = ν
            popcounts[μ,2] -= 1
            popcounts[ν,2] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1, outneighbors(g[μ],k))
                if x_μ[i] == μ
                    a[i] -= β
                elseif x_μ[i] == ν
                    a[i] += β
                end
            end
        else # recovery
            x_i[k] = 3
            popcounts[μ,2] -= 1
            popcounts[μ,3] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g[μ],k))
                a[i] -= β
            end
        end
    else # recovered node migrates
        ν = rand(outneighbors(h,μ))
        x_μ[k] = ν
        popcounts[μ,3] -= 1
        popcounts[ν,3] += 1
        a[k] = rate(k, state, mpx)
    end
end

function meanfield_fun(mpx::HeterogeneousMetaplex{SIR})
    N = nv(mpx.g[1])
    M = nv(mpx.h)
    A = adjacency_matrix.(mpx.g)
    L = laplacian_matrix(mpx.h)
    β = mpx.dynamics.β
    δ = mpx.dynamics.δ
    D = mpx.D
    f! = function(dx, x, p, t)
        for i in 1:N
            dx[1,i,:] .= .-D[1]*L'*x[1,i,:]
            dx[2,i,:] .= .-D[2]*L'*x[2,i,:]
            dx[3,i,:] .= .-D[3]*L'*x[3,i,:]
        end
        for μ in 1:M
            infection = β.*x[1,:,μ].*(A[μ]*x[2,:,μ])
            recovery = δ.*x[2,:,μ]
            dx[1,:,μ] .-= infection
            dx[2,:,μ] .+= infection .- recovery
            dx[3,:,μ] .+= recovery
        end
    end
    return f!
end
