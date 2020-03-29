
struct Metaplex{T} <: AbstractEpidemicModel where T
    g::AbstractSimpleGraph
    h::AbstractSimpleGraph
    dynamics::T
end
