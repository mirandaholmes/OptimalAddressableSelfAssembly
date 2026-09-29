using DifferentialEquations, Optimization, SciMLSensitivity, OptimizationOptimJL, CSV, DataFrames, Random, DelimitedFiles

#=
THIS FILE PRODUCES OPTIMAL BOND ENERGIES AND YIELDS AT THOSE BOND ENERGIES AT A GIVEN SET OF TIMESTEPS,
OUTPUT IS TO FILES 'YieldvsTime.txt' AND 'BindingEnergyvsTime.txt', 
NOTE THAT THE OUTPUT OF BINDING ENERGIES IS IN LEXICOGRAPHICAL ORDER
=#

# ------------------------------------ I/O and Inclusions ------------------------------------ #

# Setup Files
include("Equations.jl")
constantspath = joinpath(@__DIR__,"Constants.txt")

# Output Files
Yieldpath = joinpath(@__DIR__,"YieldvsTime.txt")
Bindingpath = joinpath(@__DIR__,"BindingEnergyvsTime.txt")

#Read Necessary Parameters from 'Constants.txt'
params = readdlm(constantspath,skipstart=2)
concentration = params[4,1]::Float64
numblocks = params[7,1]::Int64
numEdges = params[2,1]::Int64
numVertices = params[3,1]::Int64
C = Vector{Float64}(params[5,1:numblocks])
dC = Vector{Float64}(params[6,1:numblocks])

# Write output to file
function writeoutput(Yieldpath, Bindingpath, concvec, deltaarray)
    #Write Header
    maxtime = tvec[end]::Float64
    deltastring = ""
    for k in 1:numEdges
        deltastring = deltastring*"-delta$k:  "
    end
    yieldheader = "density = $concentration ; Initial deltas = $deltas ; \n Initial concentrations = $C ;\
     max time = $maxtime \n final C(Target Structure):  maxtime: \n"
    bindingheader = "density = $concentration ; Initial deltas = $deltas ; \n Initial concentrations = $C ;\
     max time = $maxtime \n"*deltastring*" maxtime: \n"
    write(Yieldpath,yieldheader)
    write(Bindingpath, bindingheader)

    #Write Data
    CSV.write(Yieldpath,Tables.table(concvec),append=true)
    CSV.write(Bindingpath,Tables.table(deltaarray),append=true)
end


# ------------------------------------ Input Parameters ------------------------------------ #

tvec = [1000.0,10000.0,100000.0] #!!! Time Steps Optimized At !!!
#Random.seed!(1234)
deltas = rand(numEdges).*-14.9
tspan = (0.0,100.0) # Time span for initial solve
maxenergy = -15.0 # Boundary condition on bonding energies

#Parameters for ODE solver and optimization
adtype = Optimization.AutoForwardDiff()
maxoptiterations = 200

#define ode problem
problem = ODEProblem(f1!,C,tspan,deltas)
sol = solve(problem,save_on = false)


# ------------------------------------ Optimization ------------------------------------ #

#Loss function
function loss(newp,t)
    newproblem = remake(problem,p = newp, tspan = (0.0,t))
    sol = solve(newproblem, Tsit5(), save_on = false)
    if Sys.free_memory() / 2^30 < 15 #Possibly helps with memory issues
        GC.gc()
    end
    loss = -sol[end,end] #Final concentration of target structure
    return loss
end

function yield_vs_Tmax(tvec,problem,adtype)
    #Initialize storage for data
    deltaarray = Array{Float64}(undef,size(tvec)[1],numEdges+1)
    concvec = Array{Float64}(undef,size(tvec)[1],2)
    i = 1 #Used for indexing
    #Optimization code for each Tmax:
    for t ∈ tvec
        optf = Optimization.OptimizationFunction((x, p) -> loss(x,t), adtype) #Define optimization loss function, in this case maximizes concentration of final structure
        optprob = Optimization.OptimizationProblem(optf, deltas, lb = ones(numEdges).*maxenergy, ub = zeros(numEdges)) #define problem and boundary conditions
        result_ode = Optimization.solve(optprob, LBFGS(), maxiters = maxoptiterations)

        for j ∈ 1:numEdges+1
            if j ≤ numEdges
                deltaarray[i,j] = -result_ode.u[j]
            else
                #deltaarray[i,:] = sort(deltaarray[i,:],rev=true) #In the case of symmetry a sorted output might make more sense
                deltaarray[i,j] = t
            end
        end
        concvec[i,1] = -result_ode.objective*numVertices/concentration
        concvec[i,2] = t
        i += 1
        println("Done optimization for t = $t")
    end
    #Write results to file
    writeoutput(Yieldpath,Bindingpath, concvec, deltaarray)
end

println("Started Parameter Optimization")
yield_vs_Tmax(tvec,problem,adtype)
println("Completed Parameter Optimization")
