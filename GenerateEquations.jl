using Graphs, Combinatorics,CSV, DataFrames, DelimitedFiles


#=
THIS FILE PRODUCES REACTION-DIFFUSION EQUATIONS AND NECESSARY CONSTANTS FOR USE IN THE OPTIMIZATION OF BOND ENERGIES AND ODE SOLVING
IT WRITES FILES NAMED 'Equations.txt' and 'Constants.txt', THE FILE EXTENSION ON EQUATIONS SHOULD BE CHANGED TO 'jl'
=#

# ------------------------------------ I/O and Inclusions ------------------------------------ #
#Includes Necessary Helpers For Set-Up
include("Helpers.jl")
#Output Paths for Setup
Equationspath = joinpath(@__DIR__,"Equations.txt")
constantspath = joinpath(@__DIR__,"Constants.txt")
# Equationspath = joinpath(@__DIR__,"../ODEParameterOptimizationAndSolving/Equations.txt")
# constantspath = joinpath(@__DIR__,"../ODEParameterOptimizationAndSolving/Constants.txt")

#Necessary Setup for Output
equationsglobal = Vector{String}([])
constantsglobal = Vector{Any}(["Constants are written in the following order; adjacency matrix, #Edges, #Vertices, density, concentration vector, derivative vector, numBlocks"])


# ------------------------------------ Constants ------------------------------------ #
#Total concentration of species in system
concentration = 0.05
#Adjacency matrix for graph
adj_mat = [0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
           1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0;
           0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0;
           0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0;
           1 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0;
           0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0;
           0 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0;
           0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 0;
           0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0;
           0 0 0 0 0 1 0 0 1 0 1 0 0 1 0 0;
           0 0 0 0 0 0 1 0 0 1 0 1 0 0 1 0;
           0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 1;
           0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0;
           0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0;
           0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1;
           0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0]
#the function in output will be named 'f<fnnum>!'
fnnum = 1

# ------------------------------------ Necessary Functions to Inititalize Setup ------------------------------------ #
#Sets indexing for deltas, should instantiate map for all adjacencies in Lexicographical order
function setDeltasmap!(deltasmap,adj_mat)
    k = 1
    for i in 1:size(adj_mat)[1]
        for j in i+1:size(adj_mat)[1]
            if adj_mat[i,j] == 1
                push!(deltasmap,[i j] => k)
                k += 1
            end
        end
    end
end


#Instantiate the list of all "Blocks", i.e. connected subgraphs of our structure
function genBlocks!(blocks,adj_mat,numVertices)
    PS = collect(powerset(1:numVertices)) #Create Vector of all possible subgraphs
    popfirst!(PS) #remove empty set
    while !isempty(PS)
        x = popfirst!(PS)
        block = minor(adj_mat,x) #Necessary to use is_connected
        if size(x)[1] == 1 || Graphs.is_connected(Graph(block)) #Check if the given block is a 'monomer' or a connected subgraph
            push!(blocks,x)
        end
    end
end


#Create dictionary mapping each block to its coordinate, and vice-versa
function genDict!(hashmap, revhashmap,blocks)
    i=0
    for x in blocks
        i += 1
        push!(hashmap,x=>i)
        push!(revhashmap,i=>x)
    end
end

#Return list of adjacencies between two blocks
function getadj(piece1,piece2,deltasmap)
    adj = []
   for i in piece1
        for j in piece2
            idx = get(deltasmap,[min(i,j) max(i,j)],0)
            if idx != 0
                push!(adj,idx)
            end
        end
    end
    return adj
end

#Function to store information about each possible interaction
#Format for Interactions is:
#In the case of a block disassembling - (sign of term in interaction, Indices of interaction energies being broken (in this case empty), Index of block dissassembling)
#In the case of a block assembling - (sign of term in interaction, Indices of interaction energies being broken, Index of first block used in assembly, Index of second block used in assembly)
function getinteractions(set,blocks,hashmap,deltasmap) #Set is the block we are generating interactions for
    vec = [] #Vector of Interactions
    #backwards interactions:
    PS = collect(powerset(set)) #Generate all possible subgraphs
    popfirst!(PS) # remove empty set
    while !isempty(PS)
        piece1 = popfirst!(PS)
        piece2 = setdiff(set,piece1)
        if piece1 in blocks && piece2 in blocks #Here we check if both blocks which assemble into set are connected !!! INCLUDE THE VALIDBOND FUNCTIONS FROM HELPERS HERE, i.e. && validbond3x3(piece1,piece2) 
            push!(vec,[1,[],get(hashmap,piece1,0),get(hashmap,piece2,0)])
            push!(vec,[-1,getadj(piece1,piece2,deltasmap),get(hashmap,set,0)])
        end
        setdiff!(PS,[piece2]) #Remove the compliment piece so as not to double count interactions
    end
    #Forwards interactions:
    setblocks = filter(b -> issubset(set,b),blocks) #setblocks is the list of blocks containing set
    while !isempty(setblocks)
        supset = popfirst!(setblocks)
        compliment = setdiff(supset,set) #compliment is the set whos disjoint union with set is supset
        if compliment in blocks #Here we check if sets compliment in supset is connected !!! INCLUDE THE VALIDBOND FUNCTIONS FROM HELPERS HERE, i.e. && validbond3x3(set,compliment)
            push!(vec,[-1,[],get(hashmap,set,0),get(hashmap,compliment,0)])
            push!(vec,[1,getadj(set,compliment,deltasmap),get(hashmap,supset,0)])
        end
    end
    return(vec)
end

#Instantiate map for each block to its list of interactions
function setInteractions!(interactionsmap,blocks,hashmap,deltasmap)
    for b in blocks
        push!(interactionsmap, b => getinteractions(b,blocks,hashmap,deltasmap))
    end
end

#  ------------------------------------ Write setup to file   ------------------------------------ #
# This code writes the ODEs generated by the interactions written above to a textfile (Equations.txt) in Julia Format
# The intention is to simply change the file extension to Equations.jl
# Also writes the necessary constants file for the solver and optimization

#add necessary constants to vector
function writeconstants!(constantsglobal,adj_mat,numEdges,numVertices,concentration,C,dC,size)
    towrite = ["",adj_mat,numEdges,numVertices,concentration,C,dC,size]
    append!(constantsglobal,towrite)
end

#Compute the vector representing all interactions
function writeinteractions!(num,equationsglobal,blocks,hashmap,deltasmap)
    equations = Vector{String}(["function f$(num)!(dC,C,δ,t)"])
    savedinteractions = Vector{String}([])
    #Loop 1 -- Precompute exponentials
    for b in blocks
        interactions = getinteractions(b,blocks,hashmap,deltasmap)
        equation = ""
        for x in interactions
            if x[2] != []
                prefix = "e_"
                suffix = ""
                for y in x[2]
                    prefix = prefix*"$y"*"_"
                    suffix = suffix*"δ[$y] + "
                end
                suffix = rstrip(suffix, [' ','+',' '])
                prefix = rstrip(prefix, ['_'])
                equation = prefix*" = exp($suffix)"
                if !(prefix in savedinteractions)
                    push!(savedinteractions,prefix)
                    push!(equations,equation)
                end
            end
        end
    end

    #Loop 2 -- Write equations
    #Here we check if the energy is one, then branch to check size x[3], otherwise branch again but use e_i_j
    for b in blocks
        idx = get(hashmap,b,0)
        equation = "dC[$idx] = "
        interactions = getinteractions(b,blocks,hashmap,deltasmap)
        for x in interactions
            if x[2] == [] #Case for interaction with energy = k_0
                if size(x)[1] == 3 #Subcase of 1 reactant
                    sgn = x[1]
                    conc = Int64(x[3])
                    if sgn == 1
                        equation = equation*"C[$conc] + "
                    else
                        equation = equation*"-1.0*C[$conc] + "
                    end
                end
                if size(x)[1] == 4 #Subcase of 2 reactants
                    sgn = x[1]
                    conc1 = Int64(x[3])
                    conc2 = Int64(x[4])
                    if sgn == 1
                        equation = equation*"C[$conc1]*C[$conc2] + "
                    else
                        equation = equation*"-1.0*C[$conc1]*C[$conc2] + "
                    end
                end
            else #Case for interaction energy k_0*exp(interactions)
                energy = "e"
                for y in x[2]
                    energy = energy*"_$y"
                end
                if size(x)[1] == 3 #Subcase for 1 reactant
                    sgn = x[1]
                    conc = Int64(x[3])
                    if sgn == 1
                        equation = equation*"C[$conc]*$energy + "
                    else
                        equation = equation*"-1.0*C[$conc]*$energy + "
                    end
                end
                if size(x)[1] == 4 #Subcase for 2 reactants
                    sgn = x[1]
                    conc1 = Int64(x[3])
                    conc2 = Int64(x[4])
                    if sgn == 1
                        equation = equation*"C[$conc1]*C[$conc2]*$energy + "
                    else
                        equation = equation*"-1.0*C[$conc1]*C[$conc2]*$energy + "
                    end
                end
            end
        end
        equation = rstrip(equation, [' ','+',' ']) #remove trailing addition signs
        push!(equations,equation)
    end
    push!(equations,"end")
    push!(equations,"")
    append!(equationsglobal,equations)
end

#Call to functions that write output
function writeoutput!(n,equationsglobal,constantsglobal,blocks,hashmap,deltasmap,adj_mat,numEdges,numVertices,concentration,C,dC)
    writeinteractions!(n,equationsglobal,blocks,hashmap,deltasmap)
    writeconstants!(constantsglobal,adj_mat,numEdges,numVertices,concentration,C,dC,size(blocks)[1])
end

#Call all necessary sub-parts of program
function setup(adj_mat)
    graph = Graph(adj_mat)
    numEdges = ne(graph)
    numVertices = nv(graph)
    blocks = []
    hashmap = Dict()
    revhashmap = Dict()
    interactionsmap = Dict()
    deltasmap = Dict()
    setDeltasmap!(deltasmap,adj_mat)
    genBlocks!(blocks,adj_mat,numVertices)
    C = zeros(Float64,size(blocks)[1])
    dC = zeros(Float64,size(blocks)[1])
    C[1:numVertices] .= concentration/numVertices
    genDict!(hashmap,revhashmap,blocks)
    setInteractions!(interactionsmap,blocks,hashmap,deltasmap)
    writeoutput!(fnnum,equationsglobal,constantsglobal,blocks,hashmap,deltasmap,adj_mat,numEdges,numVertices,concentration,C,dC)
end
setup(adj_mat) #Call the program

#Write to file
writedlm(Equationspath, equationsglobal)
writedlm(constantspath, constantsglobal)