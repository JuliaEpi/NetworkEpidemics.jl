
"""
    AbstractCompartimentalModel

An abstract type representing an arbitrary compartimental model.
"""
abstract type  AbstractCompartimentalModel end

num_states(::AbstractCompartimentalModel) = _NI("num_states")

# Metapopulation

rate(i, state, em::Metapopulation{AbstractCompartimentalModel}) = _NI("rate")

function init_rates!(a, state, mp::Metapopulation{T}) where T <: AbstractCompartimentalModel
    N = nv(mp.h)
    for k in 1:N
        a[k] = rate(k, state, mp)
    end
end

output_elem_size(mp::Metapopulation{T}) where T <: AbstractCompartimentalModel = (nv(mp.h), num_states(mp.dynamics))

function init_output(mp::Metapopulation{T}, state, nmax) where T <: AbstractCompartimentalModel
    M = nv(mp.h)
    m = num_states(mp.dynamics)
    output = [Array{Int,2}(undef, M, m) for i in 1:nmax]
    output[1] .= state
    return output
end

function update_output!(output, state, n, k, mp::Metapopulation{T}) where T <: AbstractCompartimentalModel
    output[n] .= state
end

function finalize_output(mp::Metapopulation{T}, output) where T <: AbstractCompartimentalModel
    N = nv(mp.h)
    n = length(output)
    m = num_states(mp.dynamics)
    out = [Array{Int,2}(undef, n, N) for i in 1:m]
    for i in 1:n, j in 1:m
        out[j][i,:] .= output[i][:,j]
    end
    return out
end

function finalize_meanfield(mp::Metapopulation{T}, sol) where T <: AbstractCompartimentalModel
    m = num_states(mp.dynamics)
    n = length(sol.u)
    M = nv(mp.h)
    output = [Array{Float64,2}(undef, n, M) for i in 1:m]
    for i in 1:n, j in 1:m
        output[j][i,:] .= sol.u[i][:,j]
    end
    return sol.t, output
end

# Contact Process

function init_rates!(a, state, cp::ContactProcess{T}) where T <: AbstractCompartimentalModel
    N = nv(cp.g)
    for k in 1:N
        a[k] = rate(k, state, cp)
    end
end

function init_state_mf(cp::ContactProcess{T}, x0) where T <: AbstractCompartimentalModel
    N = nv(cp.g)
    m = num_states(cp.dynamics)
    state = zeros(N,m)
    for i in 1:N
        state[i, x0[i]] = 1.0
    end
    return state
end

function finalize_meanfield(cp::ContactProcess{T}, sol) where T <: AbstractCompartimentalModel
    N = nv(cp.g)
    m = num_states(cp.dynamics)
    n = length(sol.u)
    output = [Array{Float64,2}(undef, n, N) for i in 1:m]
    for i in 1:n, j in 1:m
        output[j][i,:] .= sol.u[i][:,j]
    end
    return sol.t, output
end

# Metaplex

function init_state(mpx::Metaplex{T}, state_0) where T <: AbstractCompartimentalModel
    N = nv(mpx.g)
    M = nv(mpx.h)
    m = num_states(mpx.dynamics)
    x_i = copy(state_0[1]) # epidemic state
    x_μ = copy(state_0[2]) # location
    popcounts = zeros(Int, M, m)
    for k in 1:N
        popcounts[x_μ[k], x_i[k]] += 1
    end
    return Dict(:state_i =>x_i, :state_μ => x_μ, :popcounts => popcounts)
end

function init_rates!(a, state, mpx::Metaplex{T}) where T <: AbstractCompartimentalModel
    N = nv(mpx.g)
    for k in 1:N
        a[k] = rate(k, state, mpx)
    end
end

output_elem_size(mp::Metaplex{T}) where T <: AbstractCompartimentalModel = (nv(mp.h), num_states(mp.dynamics))

function init_output(mpx::Metaplex{T}, state, nmax) where T <: AbstractCompartimentalModel
    M = nv(mpx.h)
    m = num_states(mpx.dynamics)
    output = [Array{Int,2}(undef, M, m) for i in 1:nmax]
    output[1] .= state[:popcounts]
    return output
end

function update_output!(output, state, n, k, mpx::Metaplex{T}) where T <: AbstractCompartimentalModel
    output[n] .= state[:popcounts]
end

finalize_output(mpx::Metaplex{T}, output) where T <: AbstractCompartimentalModel = finalize_output(Metapopulation(mpx), output)

function init_state_mf(mpx::Metaplex{T}, x0) where T <: AbstractCompartimentalModel
    N = nv(mpx.g)
    M = nv(mpx.h)
    m = num_states(mpx.dynamics)
    x_i = x0[1]
    x_μ = x0[2]
    state = zeros(m,N,M)
    for i in 1:N
        state[x_i[i],i,x_μ[i]] = 1.0
    end
    return state
end

function finalize_meanfield(mpx::Metaplex{T}, sol) where T <: AbstractCompartimentalModel
    N = nv(mpx.g)
    M = nv(mpx.h)
    m = num_states(mpx.dynamics)
    n = length(sol.u)
    output = [Array{Float64,2}(undef, n, M) for i in 1:m]
    for i in 1:n, j in 1:m, k in 1:M
        output[j][i,k] = sum(sol.u[i][j,:,k])
    end
    return sol.t, output
end
