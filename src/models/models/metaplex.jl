
struct Metaplex{T} <: AbstractEpidemicModel
    g::AbstractSimpleGraph
    h::AbstractSimpleGraph
    dynamics::T
end
