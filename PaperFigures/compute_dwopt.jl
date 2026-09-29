#= 
    compute optimal weak bond strength, for given ds and T

    created Oct 31 2025

    TO DO:
    - 

=#

using DifferentialEquations, Optimization, SciMLSensitivity, OptimizationOptimJL, DelimitedFiles


whichcase = parse(Int64,ARGS[1])
ds = -parse(Float64,ARGS[2])

#whichcase = 4
#ds = -15
tlist = 10.0.^collect(3:0.25:5); #10.0.^collect(3:0.25:5);
dens = 0.05

writefile = "data_dwE"*string(whichcase)*".txt"


println("Running dwopt: whichcase = ", whichcase, ", ds = ",ds)



# optimization parameters
odesolver = AutoTsit5(Rosenbrock23())  #Tsit5()
optsolver = LBFGS();
adtype = Optimization.AutoForwardDiff()
maxoptiterations = 400

# load data specific to whichcase
if whichcase == 0
    del1 = [1,1,1,0];
    del2 = [0,0,0,1];
    nv = 4;
    nblocks = 13;
    numEdges = 4;
    include("E0.jl")
    f! = E0!
elseif whichcase==1
    del1 = [0,1,1,1,0,1,0,1,1,0,1,1]
    del2 = [1,0,0,0,1,0,1,0,0,1,0,0]
    nv = 9
    nblocks = 218
    numEdges = 12
    include("E1.jl")
    f! = E1!
elseif whichcase==2
    nv = 7
    nblocks = 95
    numEdges = 12
    del1=[0,1,1,1,1,1,0,1,0,0,0,0]
    del2=[1,0,0,0,0,0,1,0,1,1,1,1]
    include("E2.jl")
    f! = E2!
elseif whichcase==3
    nv = 8
    nblocks = 167
    numEdges = 12
    del1=[0,1,0,1,0,1,1,0,0,1,1,1]
    del2=[1,0,1,0,1,0,0,1,1,0,0,0]
    include("E3.jl")
    f! = E3!
elseif whichcase ==4
    nv = 6
    nblocks = 63
    numEdges = 15
    del1=[1,0,1,1,0,0,0,0,0,1,0,0,0,1,0]
    del2=[0,1,0,0,1,1,1,1,1,0,1,1,1,0,1]
    #del1 = [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
    #del2 = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1]
    include("E4.jl")
    f! = E4!
end


function writedeltas(dw)
    return del1*ds + del2*dw;
end



# set up ODE problem
tspan = (0.0,100.0) # Time span for initial solve 
c0 = zeros(nblocks);
c0[1:nv] .= dens/nv;   # initial concentration
problem = ODEProblem(f!,c0,tspan,writedeltas(0.1))
@time sol=solve(problem,odesolver,save_on = false)  #(forces it to compile)

#Loss function
function loss(newp,t)
    newproblem = remake(problem,p = writedeltas(newp[1]), tspan = (0.0,t))
    sol = solve(newproblem, odesolver, save_on = false)
    if Sys.free_memory() / 2^30 < 15 #Possibly helps with memory issues
        GC.gc()
    end
    loss = -sol[end,end]  #Final concentration of target structure
    return loss
end
@time loss([-1.0],100.0)  # compile, by running



# open datafile for writing (appending)
io = open(writefile, "a") 
dwopt = [];# vector for saving data


# loop through times 
for (it,tmax) in enumerate(tlist)

    println("tmax = ",tmax)

    # set up optimization problem  
    optf = Optimization.OptimizationFunction((x, p) -> loss(x,tmax), adtype) #Define optimization loss function, in this case maximizes concentration of final structure
    p0 = it==1 ?  [-0.1*rand()] : [-dwopt[end] - 0.1*rand()]  # initial dw
    optprob = Optimization.OptimizationProblem(optf, p0, lb = [ds], ub = [0]) #define problem and boundary conditions
    
    # optimize
    result_opt = Optimization.solve(optprob, optsolver, maxiters = maxoptiterations)
    
    # save data
    push!(dwopt,-result_opt.u[1])

    # print stats to screen
    println("  opt = ",-result_opt.u[1], ", niters = ",result_opt.stats.iterations)
    println("  objective = ",-result_opt.objective*nv/dens)
    println("  residual = ",result_opt.original.g_residual)
    
end

# write data to file
#writedlm(io, [-ds tlist'])
writedlm(io, [-ds dwopt'])


# close file
close(io);