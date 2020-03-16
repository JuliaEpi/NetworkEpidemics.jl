
"""
    rand_categorical(x::Vector{<:Real})

Pick one index in `1:length(x)` at random using the elements of `x` as weights.
"""
function rand_categorical(x::Vector{<:Real}, s = sum(x))
    if isempty(x)
        return 0
    end
    r = rand()
    s = sum(x)
    i = 1
    t = x[i]
    while r*s > t
        i += 1
        t += x[i]
    end
    return i
end

"""
    rand_combination(n,k)

Return a random combination of `k` elements in `1:n`.
"""
function rand_combination(n::Integer, k::Integer)
    if k < n && k > 0
        c = Vector{Int64}(undef, k)
        s = BitSet(1:n)
        for i in 1:k
            c[i] = rand(s)
            delete!(s, c[i])
        end
        return c
    else
        return []
    end
end
