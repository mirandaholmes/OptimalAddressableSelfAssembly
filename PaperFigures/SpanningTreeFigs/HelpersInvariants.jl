
include("HelpersSubTrees.jl")


function getInvariants(wc)  
    deltafile = "dataTreesE"*string(wc)*".txt"
    deltas = readdlm(deltafile,skipstart=3)
    NT = [];
    WT = [];
    for i=1:size(deltas)[1]
        NTtmp = getNT(deltas[i,:],-15.0,wc)
        push!(NT,NTtmp)
        WTtmp = WienerIndex(deltas[i,:],-15.0,wc)
        push!(WT,WTtmp)
    end
    return NT,WT
end

# I think this computes N(T) for a given set of deltas with strong interactions as deltastrong
function getNT(delta,deltastrong,whichcase)
    nlist = [9,7,8,6]
    n = nlist[whichcase]
    DFSCount(Array(2:n),ones(Int,n),zeros(Int,n),matrixify(bondenergiestoedgelistshelper(delta,deltastrong),whichcase),[[] for _ in 1:n],n)
end

# ***HELP! what do arguments represent?
function DFSCount(unvisited,InclusionTrees,ExclusionTrees,treeadjmat,children,n)  # MHC added argument "n"
    unvisited, InclusionTrees, ExclusionTrees = DFSCountHelper(1,unvisited,InclusionTrees,ExclusionTrees,treeadjmat,children,n)
    return InclusionTrees[1] + ExclusionTrees[1]
end


# ***HELP!!
function DFSCountHelper(vertex,unvisited,InclusionTrees,ExclusionTrees,treeadjmat,children,n)   # MHC added argument "n"
    for i in 1:n
        if treeadjmat[vertex,i] == 1 && i in unvisited
             filter!(e->e!=i,unvisited)
             push!(children[vertex],i)
             unvisited, InclusionTrees, ExclusionTrees = DFSCountHelper(i,unvisited,InclusionTrees,ExclusionTrees,treeadjmat,children,n)
        end
    end
    for w in children[vertex]
        InclusionTrees[vertex] *= (InclusionTrees[w] + 1)
        ExclusionTrees[vertex] += (InclusionTrees[w] + ExclusionTrees[w])
    end
    return unvisited, InclusionTrees, ExclusionTrees
end


#= 
    This is adapted from ComputeWeinerIndex.jl
=#

# MHC added for easier processing of lists of deltas
function WienerIndex(delta,deltastrong,whichcase)
    nlist = [9,7,8,6]
    n = nlist[whichcase]
    adjmat = matrixify(findall(delta.==deltastrong),whichcase)
    return WienerIndex(adjmat)
end

function WienerIndex(adj)
    dist = distance(adj)
    n = size(adj)[1]
    wi = 0;
    for i=1:n
        for j=i+1:n
            wi += dist[i,j]
        end
    end
    return wi
end

function distance(adj)
    n = size(adj)[1];
    dist = adj;
    for k=2:n
        a = adj^k;
        for i=1:n
            for j=i+1:n
                if a[i,j] >0 && dist[i,j] == 0
                    dist[i,j] = k
                    dist[j,i] = k
                end
            end
        end
    end
    return dist
end


# #Note the algorithm specified in this file is not efficient for computing the Wiener index, however, the index is only computed in small cases.
# function WienerIndex(adjmat,n)
#     wi = 0
#     for i in 1:n
#         visited = zeros(n)
#         visited[i] = 1
#         wi += WienerIndexHelper(adjmat,n,visited,[],i,1,0)
#     end
#     return wi/2
# end

# # ***HELP -- Tighe does this look correct? I changed a couple of things
# function WienerIndexHelper(adjmat,n, visited, queue, i,dist,val)
#     for j in 1:n
#         if adjmat[i,j] == 1 && visited[j] == 0   # changed from adjmat[i][j]
#             #queue.push(j)      # PROBLEM HERE: "ERROR: type Array has no field push"
#             push!(queue,j)
#             val += dist
#             visited[j] = 1
#         end
#     end
#     if isempty(queue)
#         return val
#     end
#     WienerIndexHelper(adjmat,n,visited,queue[2:end],queue[1],dist + 1,val)
# end