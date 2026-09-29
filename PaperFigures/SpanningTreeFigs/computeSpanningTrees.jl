#=
  MHC modified this from IsomorphismsClassesSubgraphs.jl
=#

using Graphs, Graphs.Experimental, DelimitedFiles 

whichcase = 4   # which example

# parameters
writepath = joinpath("dataTreesE"*string(whichcase)*".txt")
deltastrong = -15.0
deltaweak = -2.0
include("HelpersSubTrees.jl")
include("loadexample.jl")

#codedir = "/Users/mcholmes-cerfon/Library/CloudStorage/Dropbox/Work/TigheSushrut/TigheCode/PaperFigures/"


# load parameters specific to example we are looking at
n,numEdges,nblocks,adjvec = loadexample(whichcase)

# form adjacency matrix as a graph
overgraph = graph_of_advec(adjvec,n)


# compute spanning trees and write to file
function writetreesbyisoclass(whichcase::Int64)
    header = ["Data is arranged as: First row contains the number of isomorphism classes, deltaweak and deltastrong, Second row contains the number of trees per each isomorphism class, each subsequent row is the delta vectors for one of the isomorphism classes"]

    edgesets = alledgesets(whichcase)
    isotrees = isospanningtrees(whichcase)
    bondslist = edgelisttobondenergies(edgesets,deltastrong,deltaweak,numEdges)
    sortedtrees = []
    for isoclass in isotrees
        temp = bondslist
        # this seems to sort by isomorhpisms of the underlying strong bond network
        if whichcase !=4  # keep all isomorphisms for E1,E2,E3 
                          # could use all_ismorph to get isomorphisms of underlying graph, then check if these preserve edge weights
                          # not sure it's worth it though as we'll still get copies due to labelling. 
            push!(sortedtrees,filter(x -> has_isomorph(Graph(matrixify(bondenergiestoedgelistshelper(x,deltastrong),whichcase)),isoclass),temp))
        end
        if whichcase == 4  # only keep 1 copy of each isomorphism class for E4, because they're all equivalent (even with weak bonds)
            tmptrees=filter(x -> has_isomorph(Graph(matrixify(bondenergiestoedgelistshelper(x,deltastrong),whichcase)),isoclass),temp)
            push!(sortedtrees,tmptrees[1]')
        end
    end
    vals1 = [size(isotrees)[1],deltaweak,deltastrong]
    vals2 = []
    for l in sortedtrees
        push!(vals2,size(l)[1])
    end
    writedlm(writepath,header)
    open(writepath, "a") do io
        writedlm(io, [vals1])
        writedlm(io, [vals2])
        for l in sortedtrees
            writedlm(io, l)
        end
    end
end


writetreesbyisoclass(whichcase)
# bondenergies = writebondenergies(isotrees)
# println(bondenergies)