
abstract type EpidemicEvent end

struct EmptyEpidemicEvent <: EpidemicEvent

end

struct InfectionEvent <: EpidemicEvent
    node::Int32
end

struct RecoveryEvent <: EpidemicEvent
    node::Int32
end

struct CPSimpleInfectionEvent <: EpidemicEvent
    infecting_node::Int
    infected_node::Int
end

struct CPSimpleRecoveryEvent <: EpidemicEvent
    recovering_node::Int
end

struct CPSimpleRemovalEvent <: EpidemicEvent
    removed_node::Int
end
