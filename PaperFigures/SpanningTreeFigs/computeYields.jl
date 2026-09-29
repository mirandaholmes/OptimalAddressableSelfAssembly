#=
    compute yields for all the spanning trees

    created Oct 28 2025


    format of saved data in datafiles:
        [ST#, ds, dw, yields(times)]
        ... (1 row per ST) ..

=#

using DelimitedFiles, DifferentialEquations
include("loadexample.jl")


# parameters
whichcase = 4
ds = -15.0
dwlist = -[0.0,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5]
times = 10.0.^collect(3:0.25:5);
dens = 0.05;   # density
datafile = "dataYieldsE"*string(whichcase)*".txt"   # save data here

# bad pathways data
badpathwayspath = joinpath(@__DIR__,"badpathways/AllBadPathwaysE"*string(whichcase)*".txt")
badvolumespath = joinpath(@__DIR__,"badpathways/AllBadVolumesE"*string(whichcase)*".txt")
badpathways = readdlm(badpathwayspath)
badvolumes = readdlm(badvolumespath)


# open datafile for writing / appending data (write = "w"; erases previous data. append = "a")
io = open(datafile, "w") 


# spanning tree data
deltafile = joinpath(@__DIR__,"dataTreesE"*string(whichcase)*".txt")  # contains deltas of all spanning trees
ds0 = -15.0;  # ds in deltafile
dw0 = -2.0    # dw in deltafile

# read in spanning trees
deltas = readdlm(deltafile,skipstart=3)
nst = size(deltas)[1];   # total number of spanning trees


# load data specific to this case
n,numEdges,nblocks,adjvec=loadexample(whichcase)

if whichcase==1
    include("E1path.jl")
    f! = E1!
elseif whichcase==2
    include("E2path.jl")
    f! = E2!
elseif whichcase==3
    include("E3path.jl")
    f! = E3!
elseif whichcase ==4
    include("E4path.jl")
    f! = E4!
end


# compute non-pathway volume from indices and volume vector
function nonpathwayvolume(vec,badindices,badvolumes)
    val = 0
    for i in 1:size(badindices)[1]
        val += vec[badindices[i]]*badvolumes[i]/dens
    end
    return val
end



# set up ode problem, and solve it (to compile)
c0 = zeros(nblocks); # initial concentration
c0[1:n] .= dens/n;
prob = ODEProblem(f!,c0,(0,100),deltas[1,:])
odesolver = AutoTsit5(Rosenbrock23())
@time soltmp = solve(prob,odesolver) 

# define function to allow varying delta-parameters. 
# save output at t=times. Return yield vector.
function sol_params(del,tlist)
    newproblem = remake(prob,p = del, tspan = (0.0,maximum(tlist)))
    sol = solve(newproblem, odesolver, saveat=tlist)
    yield = sol[end,:] #Final concentration of target structure, at saved times
    return yield/(dens/n), sol
end


# loop through weak bonds
for dw in dwlist

    println("dw = ",dw)

    # loop through spanning trees, compute yields, and save to file
    @time for i=1:nst
        # set up delta-values
        del = deltas[i,:]
        i_ds = findall(del.==ds0)
        i_dw = findall(del.==dw0)
        del[i_ds] .= ds
        del[i_dw] .= dw

        # compute yields
        yields,sol = sol_params(del,times)

        # compute off-pathway volume
        badvol = zeros(size(times));
        for j=1:length(times)
            badvol[j] = nonpathwayvolume(sol[:,j],filter(x -> x != "",badpathways[i,:]),badvolumes[i,:]) 
        end

        # write data to file
        writedlm(io, [i ds dw yields' badvol'])
    end
end

# close file
close(io);


