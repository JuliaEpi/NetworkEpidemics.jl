
"""
    vector_array_flatten(x::Vector{Array})

Convert a vector of `N`-dimensional arrays to an `N+1`-dimensional array by concatenating along a new leading dimension.
"""
function vector_array_flatten(x::Vector{Array{T,N}}) where {T,N}
    cat(itemize.(x)..., dims=1)
end

"""
    itemize(x::Array{T,N})
    
Convert an `N`-dimensional array to an `N+1`-dimension array of the same size with a singleton leading dimension.
"""
function itemize(x::Array{T,N}) where {T,N}
    y = Array{T,N+1}(undef, 1, size(x)...)
    for i in 1:length(x)
        y[i] = x[i]
    end
    return y
end

"""
    vector_array_transpose(x, dim=1)


"""
function vector_array_transpose(x::Vector{Array{T,N}}, dim=1) where {T,N}

end
