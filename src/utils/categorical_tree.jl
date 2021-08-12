

"""
    CategoricalTree{T}
    where T <: Real

A Datastructure designed for efficiently sampling and mutating a categorical distribution with many categories.

Given an array of values, we store these, and their sums in a tree where the leaf nodes are the values, and each parent node contains the sum of up to two children, all the way up to the root node which contains the sum of all values.

This allows for ``O(1)`` read, ``O(\\log n)`` write, ``O(1)`` sum, and ``O(\\log n)`` sampling of a random index (when the values are describe a categorical distribution).
"""
struct CategoricalTree{T <: Real}
    ns::Vector{Int} # number of nodes at each layer of the tree (`ns[end]` should be 1, i.e. the root of the tree
    as::Vector{Vector{T}} # store each layer into a separate array. Each parent node contains the sum of its children
end

function CategoricalTree(a::Vector{T}) where T<:Real
    n = length(a)
    ns::Vector{Int64} = [n]
    while n > 1
        if isodd(n)
            n = Int((n+1)/2)
            push!(ns,n)
        else
            n = Int(n/2)
            push!(ns,n)
        end
    end
    as = [zeros(T,k) for k in ns]
    as[1] = a
    for i in 1:(length(ns)-1)
        for j in 2:2:ns[i]
            as[i+1][Int(j/2)] = as[i][j-1] + as[i][j]
        end
        if isodd(ns[i])
            as[i+1][ns[i+1]] = as[i][ns[i]]
        end
    end
    return CategoricalTree{T}(ns,as)
end

Base.iterate(ct::CategoricalTree) = iterate(ct.as[1])
Base.iterate(ct::CategoricalTree, state) = iterate(ct.as[1], state)

Base.isempty(ct::CategoricalTree) = isempty(ct.as[1])
Base.length(ct::CategoricalTree) = length(ct.as[1])

Base.firstindex(::CategoricalTree) = 1
Base.lastindex(ct::CategoricalTree) = length(ct)

Base.in(ct::CategoricalTree) = in(ct.as[1])

Base.getindex(ct::CategoricalTree, i) = ct.as[1][i]

function Base.setindex!(ct::CategoricalTree{T}, x::T, i) where {T}
    ct.as[1][i] = x
    j = i
    for k in 2:length(ct.ns)
        j = isodd(j) ? Int((j+1)/2) : Int(j/2)
        if isodd(ct.ns[k-1]) && j == ct.ns[k]
                ct.as[k][j] = ct.as[k-1][2*j-1]
        else
            ct.as[k][j] = ct.as[k-1][2*j-1] + ct.as[k-1][2*j]
        end
    end
end


# since we already computed the sum, we might as well make use of it
Base.sum(ct::CategoricalTree) = ct.as[end][1]

"""
    rand_categorical(ct::CategoricalTree{<:Real}, a0)

Pick one index according to a categorical distribution described by `ct`.

The leaf values `a[k]` are assumed to all be non-negative. The index `k` is picked with probability `a[k]/sum(ct)`.
"""
function rand_categorical(ct::CategoricalTree, a0 = sum(ct))
    if length(ct.ns) == 1
        return 1
    end
    r = rand()*sum(ct) # ignore a0
    i = length(ct.ns) - 1
    j = 1
    while i > 1
        if r <= ct.as[i][j]
            j = 2*j - 1
        else
            r -= ct.as[i][j]
            j = 2*j + 1
        end
        i -= 1
    end
    if r <= ct.as[i][j]
        return j
    else
        return j+1
    end
end
