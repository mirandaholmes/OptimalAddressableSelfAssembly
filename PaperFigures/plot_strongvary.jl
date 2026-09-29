# Generate data, save to file, and plot. 
# Option to generate data, or to plot. 

# created Oct 17 2025, from makedata_weakvary_diffusion.jl

# https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/

# format of datafile: 
# [whichcase dw case yielddata]
# - whichcase = example
# - tmax is time ODE solver was run for
# - case = 0 for standard model
#        = 1 for diffusion-optimal
#        = 2 for equal energies
#


using DifferentialEquations, DelimitedFiles, Plots, LaTeXStrings, ColorSchemes

include("E1.jl")
include("E2.jl")
include("E3.jl")
include("E4.jl")


# set to 0 to plot; 1 to generate data
ifgetdata = 0

# parameters
tmax = 10000.0           # time to simulate for 
dwlist = [0,1,2,4]      # weak bonds {0,1,2,3,4}
dslist = collect(9:1:20)  # strong bonds
dens = 0.05;   # density


datafile = "data_strongvary.txt"


# get colour palettes
copt1 = palette(:PuBu_6,rev=true)
copt2 = palette(:OrRd_4,rev=false)
ceq = palette(:Greens_4)


# Plot data
if ifgetdata == 0

    # set up plots
    default(fontfamily = "Computer Modern")
    default(linewidth=2,titlefontsize=12,legendfontsize=8,legendtitlefontsize =10)
    #default(linewidth=2,titlefontsize=14,legendfontsize=10,legendtitlefontsize =12)

    # set up subplots
    l = @layout grid(2,2)  #[grid(2,2) a{0.2w}]
    plt = plot(layout=l) 

    # add title/axis labels
    for i in 1:4
        plot!(xaxis=L"\delta_s",yaxis="yield",ylimits=(0,1),title="E"*string(i)*", "*L"T=10^4",subplot=i)   
    end

    # # set up legend
    # plot!(plt, showaxis = false, grid = false, leg_title=L"\delta_w", legend=(-0.4,0.7),  subplot=5)
    # for it=1:3
    #     plot!(plt,[0; 100],[NaN; NaN],label=string(dwlist[it]),color=copt1[it+1], subplot=5)
    # end

    # load data
    data = readdlm(datafile)

    for (sp,whichcase) in enumerate([1,2,3,4])

        println("whichcase = ",whichcase);

        # loop through and load data and plot it
        for (iw,dw) in enumerate(dwlist)
            rowind = findall((data[:,1] .== whichcase) .& (data[:,2] .== dw) .& (data[:,3] .== 0))
            yield = data[rowind,4:end]';
            plot!(plt,dslist,yield,color=copt1[iw+1],label=string(dw),leg_title=L"\delta_w",subplot=sp)
        end
        
        # remove legends in all but one plot
        if whichcase != 1
            plot!(plt,legend=false,subplot=sp) 
        end
    end


    # display plot
    display(plt)

    # save (adjust size until font sizes are large enough when saved)
    #plot!(plt,size=(600,400))
    #savefig("./strongvary.pdf")  
end


# Get the data
if ifgetdata == 1

    # open datafile for writing
    io = open(datafile, "w") 

    for (sp,whichcase) in enumerate([1,2,3,4])

        println("whichcase = ",whichcase);

        if whichcase == 0
            del1 = [1,1,1,0];
            del2 = [0,0,0,1];
            nv = 4;
            nblocks = 13;
            f! = E0!
        elseif whichcase==1
            del1 = [0,1,1,1,0,1,0,1,1,0,1,1]
            del2 = [1,0,0,0,1,0,1,0,0,1,0,0]
            nv = 9
            nblocks = 218
            f! = E1!
        elseif whichcase==2
            nv = 7
            nblocks = 95
            del1=[0,1,1,1,1,1,0,1,0,0,0,0]
            del2=[1,0,0,0,0,0,1,0,1,1,1,1]
            f! = E2!
        elseif whichcase==3
            nv = 8
            nblocks = 167
            del1=[0,1,0,1,0,1,1,0,0,1,1,1]
            del2=[1,0,1,0,1,0,0,1,1,0,0,0]
            f! = E3!
        elseif whichcase ==4
            nv = 6
            nblocks = 63
            del1=[1,0,1,1,0,0,0,0,0,1,0,0,0,1,0]
            del2=[0,1,0,0,1,1,1,1,1,0,1,1,1,0,1]
            f! = E4!
        end

        # initial concentration
        c0 = zeros(nblocks);
        c0[1:nv] .= dens/nv;

        # set up ode problem
        prob = ODEProblem(f!,c0,(0,10),-del1)

        # solve it
        odesolver = AutoTsit5(Rosenbrock23())
        @time sol = solve(prob,odesolver)  


        # define function to allow varying parameters
        function sol_params(deltas,t)
            newproblem = remake(prob,p = -deltas, tspan = (0.0,t))
            sol = solve(newproblem, odesolver, save_on = false)
            yield = sol[end,end] #Final concentration of target structure
            return yield/(dens/nv)
        end


        # loop through different weak bond strengths
        for (iw,dw) in enumerate(dwlist)

            println("dw = ",dw)

            # Solve for yield as a function of strong bonds
            yield = zeros(size(dslist))
            @time for i in 1:length(dslist)
                ds = dslist[i]
                deltas = ds*del1 + dw*del2   # bond strengths
                yield[i] = sol_params(deltas,tmax)
            end
            # write data to file
            writedlm(io, [whichcase dw 0 yield'])

        end  # loop through tmax
    end  # loop through whichcase
    # close file
    close(io);
end