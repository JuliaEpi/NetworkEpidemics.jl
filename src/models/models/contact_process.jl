
struct ContactProcess{T} <: AbstractEpidemicModel where T
    g::AbstractSimpleGraph
    dynamics::T
end
