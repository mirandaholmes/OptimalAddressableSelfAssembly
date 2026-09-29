# MATRIX MINOR FUNCTION:
function minor(A, vec)
    n = size(vec)[1]
    B = zeros(Int64,n,n)
    for i in 1:n
        for j in 1:n
            B[i,j] = A[vec[i],vec[j]]
        end
    end
    return B
end

# FUNCTIONS TO PREVENT OBSTRUCTED BONDS:

# 3x3 CASE:
function validbond3x3(b1,b2)
    if b1 == [5]
        return !issubset([2,4,6,8],b2)
    elseif b2 == [5]
        return !issubset([2,4,6,8],b1)
    else
        return true
    end
end

# HEX GRID 7 VERTEX CASE
function validbond7(b1,b2)
    if b1 == [1]
        return b2 != [2,3,4,5,6,7]
    elseif b2 == [1]
        return b1 != [2,3,4,5,6,7]
    else
        return true
    end
end

# 4x4 CASE
#Possible inner obstructed species
innerlist4x4 = [[6,7,10,11],[6,7,10],[6,7,11],[6,10,11],[7,10,11],[6,7],[6,10],[7,11],[10,11],[6],[7],[10],[11]]
#Possible outer obstructed species
outerlist4x4 = [[2,3,5,8,9,12,14,15],[2,3,5,8,9,11,14],[2,3,5,8,10,12,15],[2,5,7,9,12,14,15],[3,6,8,9,12,14,15],[2,3,5,8,10,11],[2,5,7,9,11,14],[3,6,8,10,12,15],[6,7,9,12,14,15],
[2,5,7,10],[3,6,8,11],[6,9,11,14],[7,10,12,15]]
function validbond4(b1,b2)
    if b1 in innerlist4x4
        return isempty(filter(x -> issubset(x,b2),outerlist4x4))
    elseif b2 in innerlist4x4
        return isempty(filter(x -> issubset(x,b1),outerlist4x4))
    else
        return true
    end
end

# ONE POSSIBLY USEFUL FUNCTION I USED FOR DEBUGGING
# Output is in the form of two columns representing the (i,j) adjacencies of each interaction
#interactionspath = joinpath(@__DIR__,"../ODEParameterOptimizationAndSolving/Interactions.txt")
function checkinteractions(adj_mat)
    interactions = Vector{Tuple{Int64,Int64}}([])
    for i in 1:size(adj_mat)[1]
        for j in i+1:size(adj_mat)[1]
            if adj_mat[i,j] == 1
                push!(interactions,(i,j))
            end
        end
    end
    writedlm(interactionspath,interactions)
end