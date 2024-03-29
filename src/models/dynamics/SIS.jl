
"""
    SIS

A type representing a Susceptible-Infected-Susceptible model with infection rate `β` and recovery rate `γ`.
"""
struct SIS <: AbstractCompartimentalModel
    β::Real
    γ::Real
end

SI(sis::SIS) = SI(sis.β)

num_states(::SIS) = 2

# Metapopulation

function rate(i, state, mp::Metapopulation{SIS})
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    D = mp.D
    h = mp.h
    od = outdegree(h,i)
    β*state[i,1]*state[i,2] + γ*state[i,2] + od * sum(D .* state[i,:])
end

function update_state_and_rates!(state, a, k, mp::Metapopulation{SIS})
    N = nv(mp.h)
    D = mp.D
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    h = mp.h
    od = outdegree(h,k)
    p = vcat(β*state[k,1]*state[k,2], γ*state[k,2], D.*state[k,:]*od)
    j = sample(ProbabilityWeights(p, a[k]))
    if j == 1
        state[k,1] -= 1
        state[k,2] += 1
        a[k] = rate(k, state, mp)
    elseif j == 2
        state[k,1] += 1
        state[k,2] -= 1
        a[k] = rate(k, state, mp)
    else
        μ = k
        ν = rand(outneighbors(h, μ))
        state[μ,j-2] -= 1 # NB. j-2 is the state of the individual that migrates
        state[ν,j-2] += 1 # NB. because the first 2 elements of p are the reactions
        a[μ] = rate(μ, state, mp)
        a[ν] = rate(ν, state, mp)
    end
end

function meanfield_fun(mp::Metapopulation{SIS})
    L = laplacian_matrix(mp.h) |> float
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    D = mp.D
    f! = function(dx, x, p, t)
        tmp = (β.*x[:,1] .- γ).*x[:,2]
        dx[:,1] .= .-tmp .- D[1].*(L'*x[:,1])
        dx[:,2] .= tmp .- D[2].*(L'x[:,2])
    end
    return f!
end

# Contact Process

function rate(k, state, cp::ContactProcess{SIS})
    if state[k] == 1
        l = length(filter(i->state[i]==2, inneighbors(cp.g, k)))
        return cp.dynamics.β * l
    else
        return cp.dynamics.γ
    end
end


function update_state_and_rates!(state, a, k, cp::ContactProcess{SIS})
    g = cp.g
    β = cp.dynamics.β
    if state[k] == 1
        state[k] = 2
        for i in Iterators.filter(j->state[j]==1, outneighbors(g, k))
            a[i] += β
        end
    else
        state[k] = 1
        for i in Iterators.filter(j->state[j]==1, outneighbors(g, k))
            a[i] -= β
        end
    end
    a[k] = rate(k, state, cp)
end

function update_output!(output, state, n, k, ::ContactProcess{SIS})
    if state[k] == 2 # newly infected node
        output[n] = output[n-1] + [-1, 1]
    else
        output[n] = output[n-1] + [1, -1]
    end
end


function init_state_mf(::ContactProcess{SIS}, x0)
    return float(x0 .- 1)
end

function meanfield_fun(cp::ContactProcess{SIS})
    A = adjacency_matrix(cp.g) |> float
    β = cp.dynamics.β
    γ = cp.dynamics.γ
    f! = function(dx, x, p, t)
        dx .= β.*(1.0 .- x).*(A*x) .- γ.*x
    end
    return f!
end

function finalize_meanfield(cp::ContactProcess{SIS}, sol)
    N = nv(cp.g)
    n = length(sol.u)
    output = Array{Float64,2}(undef, n, N)
    for i in 1:n
        output[i,:] .= sol.u[i]
    end
    return sol.t, output
end

# Metaplex

function rate(k, state, mpx::Metaplex{SIS})
    β = mpx.dynamics.β
    γ = mpx.dynamics.γ
    D = mpx.D
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    od = outdegree(mpx.h, x_μ[k])
    if x_i[k] == 1
        l = length(filter(i->x_i[i]==2 && x_μ[i]==x_μ[k], inneighbors(mpx.g, k)))
        return β*l + D[1]*od
    else
        return γ + D[2]*od
    end
end

function update_state_and_rates!(state, a, k, mpx::Metaplex{SIS})
    N = nv(mpx.g)
    M = nv(mpx.h)
    g = mpx.g
    h = mpx.h
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    od = outdegree(mpx.h, x_μ[k])
    popcounts = state[:popcounts]
    β = mpx.dynamics.β
    γ = mpx.dynamics.γ
    D = mpx.D
    μ = x_μ[k]
    if x_i[k] == 1
        if rand()*a[k] < D[1]*od # susceptible node migrates
            ν = rand(outneighbors(h, μ))
            x_μ[k] = ν
            popcounts[μ,1] -= 1
            popcounts[ν,1] += 1
            a[k] = rate(k, state, mpx)
        else # susceptible node becomes infected
            x_i[k] = 2
            popcounts[μ,1] -= 1
            popcounts[μ,2] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g,k))
                a[i] += β
            end
        end
    else
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
        else # infected node recovers
            x_i[k] = 1
            popcounts[μ,2] -= 1
            popcounts[μ,1] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g,k))
                a[i] -= β
            end
        end
    end
end

function meanfield_fun(mpx::Metaplex{SIS})
    N = nv(mpx.g)
    M = nv(mpx.h)
    A = adjacency_matrix(mpx.g) |> float
    L = laplacian_matrix(mpx.h) |> float
    β = mpx.dynamics.β
    γ = mpx.dynamics.γ
    D = mpx.D
    f! = function(dx, x, p, t)
        for i in 1:N
            dx[1,i,:] .= .-D[1]*L'*x[1,i,:]
            dx[2,i,:] .= .-D[2]*L'*x[2,i,:]
        end
        for μ in 1:M
            infection = β.*x[1,:,μ].*(A*x[2,:,μ]) .- γ.*x[2,:,μ]
            dx[1,:,μ] .-= infection
            dx[2,:,μ] .+= infection
        end
    end
    return f!
end

# Heterogeneous Metaplex

function rate(k, state, mpx::HeterogeneousMetaplex{SIS})
    β = mpx.dynamics.β
    γ = mpx.dynamics.γ
    D = mpx.D
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    μ = x_μ[k]
    od = outdegree(mpx.h, x_μ[k])
    if x_i[k] == 1
        l = length(filter(i->x_i[i]==2 && x_μ[i]==x_μ[k], inneighbors(mpx.g[μ], k)))
        return β*l + D[1]*od
    else
        return γ + D[2]*od
    end
end

function update_state_and_rates!(state, a, k, mpx::HeterogeneousMetaplex{SIS})
    g = mpx.g
    h = mpx.h
    x_i = state[:state_i]
    x_μ = state[:state_μ]
    popcounts = state[:popcounts]
    β = mpx.dynamics.β
    γ = mpx.dynamics.γ
    D = mpx.D
    μ = x_μ[k]
    od = outdegree(mpx.h, x_μ[k])
    if x_i[k] == 1
        if rand()*a[k] < D[1]*od # susceptible node migrates
            ν = rand(outneighbors(h, μ))
            x_μ[k] = ν
            popcounts[μ,1] -= 1
            popcounts[ν,1] += 1
            a[k] = rate(k, state, mpx)
        else # susceptible node becomes infected
            x_i[k] = 2
            popcounts[μ,1] -= 1
            popcounts[μ,2] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g[μ],k))
                a[i] += β
            end
        end
    else
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
        else # infected node recovers
            x_i[k] = 1
            popcounts[μ,2] -= 1
            popcounts[μ,1] += 1
            a[k] = rate(k, state, mpx)
            for i in Iterators.filter(j->x_i[j]==1 && x_μ[j]==μ, outneighbors(g[μ],k))
                a[i] -= β
            end
        end
    end
end

function meanfield_fun(mpx::HeterogeneousMetaplex{SIS})
    N = nv(mpx.g[1])
    M = nv(mpx.h)
    A = mpx.adj
    L = laplacian_matrix(mpx.h) |> float
    β = [dyn.β for dyn in mpx.dynamics]
    γ = [dyn.γ for dyn in mpx.dynamics]
    D = mpx.D
    infection = Vector{Float64}(undef, N)
    f! = function(dx, x, p, t)
        for i in 1:N
            dx[1,i,:] .= .-D[1]*L'*view(x,1,i,:)
            dx[2,i,:] .= .-D[2]*L'*view(x,2,i,:)
        end
        for μ in 1:M
            infection .= β[μ].*view(x,1,:,μ).*(A[μ]*view(x,2,:,μ)) .- γ[μ].*view(x,2,:,μ)
            dx[1,:,μ] .-= infection
            dx[2,:,μ] .+= infection
        end
    end
    return f!
end

# HeterogeneousMetapopulation

function rate(i, state, mp::HeterogeneousMetapopulation{SIS})
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    D = mp.D
    h = mp.h
    ks = mp.ks
    N = mp.N
    βᵢ = β * ks[i] / N
    return βᵢ*state[i,1]*state[i,2] + γ*state[i,2] + sum(D .* state[i,:]) * (outdegree(h,i))
end

function update_state_and_rates!(state, a, k, mp::HeterogeneousMetapopulation{SIS})
    M = nv(mp.h)
    N = mp.N
    ks = mp.ks
    D = mp.D
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    h = mp.h
    βₖ = β * ks[k] / N
    p = vcat(βₖ*state[k,1]*state[k,2], γ*state[k,2], D.*state[k,:]*(outdegree(h,k)))
    j = sample(ProbabilityWeights(p, a[k]))
    if j == 1
        state[k,1] -= 1
        state[k,2] += 1
        a[k] = rate(k, state, mp)
    elseif j == 2
        state[k,1] += 1
        state[k,2] -= 1
        a[k] = rate(k, state, mp)
    else
        μ = k
        ν = rand(outneighbors(h, μ))
        state[μ,j-2] -= 1 # NB. j-2 is the state of the individual that migrates
        state[ν,j-2] += 1 # NB. because the first 2 elements of p are the reactions
        a[μ] = rate(μ, state, mp)
        a[ν] = rate(ν, state, mp)
    end
end

function meanfield_fun(mp::HeterogeneousMetapopulation{SIS})
    #L = normalized_laplacian(mp.h)
    L = laplacian_matrix(mp.h) |> float
    ks = mp.ks
    N = mp.N
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    D = mp.D
    βs = β * ks / N
    f! = function(dx, x, p, t)
        tmp = (βs .* x[:,1] .- γ).*x[:,2]
        dx[:,1] .= .-tmp .- D[1].*(L'*x[:,1])
        dx[:,2] .= tmp .- D[2].*(L'x[:,2])
    end
    return f!
end
