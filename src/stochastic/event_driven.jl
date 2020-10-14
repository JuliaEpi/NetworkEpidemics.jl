
function event_driven(cp::ContactProcess{SIR}, x₀; tmax=100.0)
    N = nv(cp.g)
    β = cp.dynamics.β
    γ = cp.dynamics.δ
    events = PriorityQueue{EpidemicEvent,Float64}()
    # maximum 2 events per node, so 2*N events max
    ts = Vector{Float64}(undef, 2*N)
    S = Vector{Int32}(undef, 2*N)
    I = Vector{Int32}(undef, 2*N)
    R = Vector{Int32}(undef, 2*N)
    k = 1
    ts[k] = 0.0
    S[k] = N
    I[k] = 0
    R[k] = 0
    state = copy(x₀)
    for u in 1:N
        if x₀[u] == 2
            events[InfectionEvent(u)] = 0.0
        end
    end
    while !isempty(events)
        (e,t) = dequeue_pair!(events)
        if t != ts[k]
            k += 1
            ts[k] = t
            S[k] = S[k-1]
            I[k] = I[k-1]
            R[k] = R[k-1]
        end
        if typeof(e) == InfectionEvent
            S[k] -= 1
            I[k] += 1
            u = e.node
            state[u] = 2
            if γ != 0
                τᵣ = log(1.0/rand())/γ
                events[RecoveryEvent(u)] = t+τᵣ
            else
                τᵣ = Inf
            end
            for v in filter(x->state[x]==1, neighbors(cp.g, u))
                τᵢ = log(1.0/rand())/β
                if τᵢ >= τᵣ
                    continue
                end
                if haskey(events, InfectionEvent(v))
                    if t+τᵢ >= events[InfectionEvent(v)]
                        continue
                    end
                end
                events[InfectionEvent(v)] = t + τᵢ
            end
        else # recovery
            I[k] -= 1
            R[k] += 1 
        end
    end
end