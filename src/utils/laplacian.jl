
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
