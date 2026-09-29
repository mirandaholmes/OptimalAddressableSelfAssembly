using DifferentialEquations, Optimization, SciMLSensitivity, OptimizationOptimJL, CSV, DataFrames, DelimitedFiles, Random

#add Optimization SciMLSensitivity OptimizationOptimJL CSV DataFrames DelimitedFiles Random
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
Bindingpath = joinpath(@__DIR__,"params1.txt")

#Read Necessary Parameters from 'Constants.txt'
params = readdlm(constantspath,skipstart=2)
concentration = params[4,1]::Float64
numblocks = params[7,1]::Int64
numEdges = params[2,1]::Int64
numVertices = params[3,1]::Int64
C = Vector{Float64}(params[5,1:numblocks])
dC = Vector{Float64}(params[6,1:numblocks])

# Write output to file
function writeoutput(Bindingpath, concvec, deltaarray,t)
    # Write Data
    # Format: density, time, c4, gradient, optimal deltas
    delinfo = [concentration;t;concvec';round.(deltaarray',digits=5)]';
    CSV.write(Bindingpath,Tables.table(delinfo),append=true)
end


# ------------------------------------ Input Parameters ------------------------------------ #

tvec = [10000.0,10000.0,10000.0]; #!!! Time Steps Optimized At !!!
maxenergy = -15.0 # Boundary condition on bonding energies

#Parameters for ODE solver and optimization
maxoptiterations = 200
odesolver = Tsit5();
optsolver = LBFGS();

#define ode problem
deltas = rand(numEdges).*maxenergy  # Initialization deltas 
tspan = (0.0,100.0) # Time span for initial solve 
problem = ODEProblem(f1!,C,tspan,deltas)
@time solve(problem,odesolver,save_on = false)  #(forces it to compile)



# ------------------------------------ Optimization Functions ------------------------------------ #

#Loss function
function loss(newp,t)
    newproblem = remake(problem,p = newp, tspan = (0.0,t))
    sol = solve(newproblem, odesolver, save_on = false)
    loss = -sol[end,end] #Final concentration of target structure
    return loss
end

@time loss(rand(numEdges).*maxenergy,100.0)  # compile, by running



function yield_vs_Tmax(tvec)
    #Initialize storage for data
    deltaarray = Array{Float64}(undef,1,numEdges)
    concvec = Array{Float64}(undef,1,2)
    #Optimization code for each Tmax:
    for t ∈ tvec
        println("Starting optimization for t = $t")
        optf = Optimization.OptimizationFunction((x, p) -> loss(x,t), Optimization.AutoForwardDiff()) #Define optimization loss function, in this case maximizes concentration of final structure
        deltas = rand(numEdges).*maxenergy  # Initialization deltas
        optprob = Optimization.OptimizationProblem(optf, deltas, lb = ones(numEdges).*maxenergy, ub = zeros(numEdges)) #define problem and boundary conditions
        @time result_ode = Optimization.solve(optprob, optsolver, maxiters = maxoptiterations)

        for j ∈ 1:numEdges
            deltaarray[j] = -result_ode.u[j]
        end
        concvec[1] = -result_ode.objective*numVertices/concentration
        concvec[2] = result_ode.original.g_residual
        println("Done optimization for t = $t")
        #Write results to file
        writeoutput(Bindingpath, concvec, deltaarray, t)
    end
end


# ------------------------------------ Optimize ------------------------------------ #

println("Started Parameter Optimization")
yield_vs_Tmax(tvec)
println("Completed Parameter Optimization")




# ------------------------------------ Debugging ------------------------------------ #

# debugging: test autodiff
#using ForwardDiff
#function f(x)
#    _prob = remake(problem, p=x[1:end-1], tspan=(0.0, x[end]))
#    solve(_prob, reltol=1e-6, abstol=1e-6, saveat=1)[4, end]
#end
#x = [deltas; 10.0]
#@time dx = ForwardDiff.gradient(f, x)
#println("dx = ",dx);


#=
# Simple function to check conservation of mass
function conservesmass(soln)
    vec = zeros(size(soln)[2])
    for i in [1:size(soln)[2]]
        for b in blocks
            vec[i] += soln[get(hashmap,b,0),i]*size(b)[1]
        end
    end
    return vec
end
result = conservesmass(sol)
println("masses: $result")
=#
