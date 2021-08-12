"""
    normalized_laplacian(g)

Return the random-walk normalized laplacian of `g`, i.e. the matrix ``D^{-1}L``, where ``D`` is the diagonal degree matrix, and ``L`` is the standard laplacian matrix.
"""
function normalized_laplacian(g::SimpleGraph)
    L = float(laplacian_matrix(g))
    for i in vertices(g)
        if L[i,i] != 0
            L[i,:] /= L[i,i]
        end
    end
    return L
end

function normalized_laplacian(g::SimpleDiGraph)
    L = float(laplacian_matrix(g))
    for i in vertices(g)
        if L[i,i] != 0
            L[i,:] /= L[i,i]
        end
    end
    return L
end
