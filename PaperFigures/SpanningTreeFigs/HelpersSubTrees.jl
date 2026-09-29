using Combinatorics, Graphs

#= 
    these functions were taken from IsomorphismsClassesSubgraphs.jl
=#

# form a graph from list of edges (vector of 1/0s)
function graph_of_advec(adjvec,n)
    adjmat = zeros(Int,n,n)
    for i in 1:n
        adjmat[i,:] = adjvec[1,n*(i-1) + 1 : n*i]
    end
    overgraph = Graph(adjmat)  # convert to graph structure
    return overgraph
end

# form an adjacency matrix from list of edges (vector of 1/0s)
function adjmat_of_advec(adjvec,n)
    adjmat = zeros(Int,n,n)
    for i in 1:n
        adjmat[i,:] = adjvec[1,n*(i-1) + 1 : n*i]
    end
    return adjmat
end

#Filter the list of trees to have at most one representative for each isomorphism class
# ***HELP!!! -- what kind of isomorphism?
function filteriso!(trees)
    i = 0
    while i < size(trees)[1]
        i += 1
        filter!(x -> x == trees[i] || !has_isomorph(trees[i],x),trees)
    end
end

#Check if given tree is included in the given target structure
function issubgraph(tree)
    for e in edges(tree)
        if !has_edge(overgraph,e)
            return false
        end
    end
    return true
end

# ***HELP
function isospanningtrees(whichcase)
    isotrees = map(x -> Graph(matrixify(x,whichcase)),alledgesets(whichcase))
    filteriso!(isotrees)
    return isotrees
end

# ***HELP
function writebondenergies(trees,overgraph,deltastrong,deltaweak)
    j = 1
    bondenergiesvec = Vector{Vector{Float64}}(undef, size(trees)[1])
    for t in trees
        deltas = Vector{Float64}(undef,ne(overgraph))
        i = 1
        for e in edges(overgraph)
            if e in edges(t)
                deltas[i] = deltastrong
            else
                deltas[i] = deltaweak
            end
            i+=1
        end
        bondenergiesvec[j] = deltas
        j+=1
    end
    return bondenergiesvec
end

# ***HELP!!!!
function bondenergiestoedgelists(bondenergies,deltastrong)
   return map(b -> bondenergiestoedgelistshelper(b,deltastrong), bondenergies)
end

# ***HELP!!!!
function bondenergiestoedgelistshelper(be,deltastrong)
   i = 1
   res = []
    for i in 1:size(be)[1]
        if be[i] == deltastrong
            push!(res,i)
        end
        i += 1
    end
    return res
end

# ***HELP!!!!
function edgelisttobondenergies(edgelist,deltastrong,deltaweak,numEdges)
    return map(e -> edgelisttobondenergieshelper(e,deltastrong,deltaweak,numEdges),edgelist)
end

# ***HELP!!!!
function edgelisttobondenergieshelper(el,deltastrong, deltaweak,numEdges)
    deltas = ones(numEdges).*deltaweak
    for idx in el
        deltas[idx] = deltastrong
    end
    return deltas
end



#=
    these functions were taken from AllSubTrees.jl
=#

# produces an adjacency matrix from an edgeset
# ***HELP!!! what is format of edgeset?
function matrixify(edgeset,whichcase::Int64)
    if whichcase == 1
        res = matrix1(edgeset)
    end
    if whichcase == 2
        res = matrix2(edgeset)
    end
    if whichcase == 3
        res = matrix3(edgeset)
    end
    if whichcase == 4
        res = matrix4(edgeset)
    end

    return res
end

# ***HELP!!!
function alledgesets(whichcase::Int64)
    sizes = [8,6,7,5]  # ***HELP -- are these sizes correct? what do they represent? is it supposed to be n-1?
    candidates = collect(powerset(1:12))
    edgesets = filter(x -> size(x)[1]==sizes[whichcase], candidates)
    return filter(e -> is_connected(Graph(matrixify(e,whichcase))),edgesets)
end



# functions specific to each example, for forming adjacency matrix
function matrix1(edgeset)
    res = zeros(9,9)
    if  1 in edgeset
        res[1,2] = 1
        res[2,1] = 1
    end
    if 2 in edgeset
        res[1,4] = 1
        res[4,1] = 1
    end
    if 3 in edgeset
        res[2,3] = 1
        res[3,2] = 1
    end
    if 4 in edgeset
        res[2,5] = 1
        res[5,2] = 1
    end
    if 5 in edgeset
        res[3,6] = 1
        res[6,3] = 1
    end
    if 6 in edgeset
        res[4,5] = 1
        res[5,4] = 1
    end
    if 7 in edgeset
        res[4,7] = 1
        res[7,4] = 1
    end
    if 8 in edgeset
        res[5,6] = 1
        res[6,5] = 1
    end
    if 9 in edgeset
        res[5,8] = 1
        res[8,5] = 1
    end
    if 10 in edgeset
        res[6,9] = 1
        res[9,6] = 1
    end
    if 11 in edgeset
        res[7,8] = 1
        res[8,7] = 1
    end
    if 12 in edgeset
        res[8,9] = 1
        res[9,8] = 1
    end
    return res
end

function matrix2(edgeset)
    res = zeros(7,7)
    if  1 in edgeset
        res[1,2] = 1
        res[2,1] = 1
    end
    if 2 in edgeset
        res[1,3] = 1
        res[3,1] = 1
    end
    if 3 in edgeset
        res[1,4] = 1
        res[4,1] = 1
    end
    if 4 in edgeset
        res[1,5] = 1
        res[5,1] = 1
    end
    if 5 in edgeset
        res[1,6] = 1
        res[6,1] = 1
    end
    if 6 in edgeset
        res[1,7] = 1
        res[7,1] = 1
    end
    if 7 in edgeset
        res[2,3] = 1
        res[3,2] = 1
    end
    if 8 in edgeset
        res[2,7] = 1
        res[7,2] = 1
    end
    if 9 in edgeset
        res[3,4] = 1
        res[4,3] = 1
    end
    if 10 in edgeset
        res[4,5] = 1
        res[5,4] = 1
    end
    if 11 in edgeset
        res[5,6] = 1
        res[6,5] = 1
    end
    if 12 in edgeset
        res[6,7] = 1
        res[7,6] = 1
    end
    return res
end

function matrix3(edgeset)
        res = zeros(8,8)
    if  1 in edgeset
        res[1,2] = 1
        res[2,1] = 1
    end
    if 2 in edgeset
        res[1,3] = 1
        res[3,1] = 1
    end
    if 3 in edgeset
        res[1,5] = 1
        res[5,1] = 1
    end
    if 4 in edgeset
        res[2,4] = 1
        res[4,2] = 1
    end
    if 5 in edgeset
        res[2,6] = 1
        res[6,2] = 1
    end
    if 6 in edgeset
        res[3,4] = 1
        res[4,3] = 1
    end
    if 7 in edgeset
        res[3,7] = 1
        res[7,3] = 1
    end
    if 8 in edgeset
        res[4,8] = 1
        res[8,4] = 1
    end
    if 9 in edgeset
        res[5,6] = 1
        res[6,5] = 1
    end
    if 10 in edgeset
        res[5,7] = 1
        res[7,5] = 1
    end
    if 11 in edgeset
        res[6,8] = 1
        res[8,6] = 1
    end
    if 12 in edgeset
        res[7,8] = 1
        res[8,7] = 1
    end
    return res
end

function matrix4(edgeset)
        res = zeros(6,6)
    if  1 in edgeset
        res[1,2] = 1
        res[2,1] = 1
    end
    if 2 in edgeset
        res[1,3] = 1
        res[3,1] = 1
    end
    if 3 in edgeset
        res[1,4] = 1
        res[4,1] = 1
    end
    if 4 in edgeset
        res[1,5] = 1
        res[5,1] = 1
    end
    if 5 in edgeset
        res[1,6] = 1
        res[6,1] = 1
    end
    if 6 in edgeset
        res[2,3] = 1
        res[3,2] = 1
    end
    if 7 in edgeset
        res[2,4] = 1
        res[4,2] = 1
    end
    if 8 in edgeset
        res[2,5] = 1
        res[5,2] = 1
    end
    if 9 in edgeset
        res[2,6] = 1
        res[6,2] = 1
    end
    if 10 in edgeset
        res[3,4] = 1
        res[4,3] = 1
    end
    if 11 in edgeset
        res[3,5] = 1
        res[5,3] = 1
    end
    if 12 in edgeset
        res[3,6] = 1
        res[6,3] = 1
    end
    if 13 in edgeset
        res[4,5] = 1
        res[5,4] = 1
    end
    if 14 in edgeset
        res[4,6] = 1
        res[6,4] = 1
    end
    if 15 in edgeset
        res[5,6] = 1
        res[6,5] = 1
    end
    return res
end