import Base.length
import Base.isempty
import Base.getindex
import Base.setindex!
import Base.firstindex
import Base.lastindex
import Base.sum

# Use a binary tree to describe a categorical distribution
# This allows for O(1) read, O(log n) write and O(log n) generation of a random index
# Optimized for categorical distributions with a large number of possible outcomes

struct CategoricalTree{T<:Real}
    ns::Vector{<:Integer} # number of nodes at each layer of the tree (`ns[end]` should be 1, i.e. the root of the tree
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

isempty(ct::CategoricalTree{<:Real}) = isempty(ct.as[1])
length(ct::CategoricalTree{<:Real}) = length(ct.as[1])

firstindex(::CategoricalTree) = 1
lastindex(ct::CategoricalTree) = length(ct)

getindex(ct::CategoricalTree{<:Real}, i) = ct.as[1][i]

function setindex!(ct::CategoricalTree{<:Real}, x, i)
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
sum(ct::CategoricalTree{<:Real}) = ct.as[end][1]

"""
    rand_categorical(ct::CategoricalTree{<:Real}, a0)

Pick one index according to a categorical distribution described by `ct`.
"""
function rand_categorical(ct::CategoricalTree{<:Real}, a0 = sum(ct))
    r = rand()*sum(ct)
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
