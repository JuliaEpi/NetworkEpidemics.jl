
struct SI{T} <: AbstractEpidemicModel
    β::Real
    model::T
end

function init_state(mp::SI{Metapopulation}, state_0, t)
    
end
