# Plot yields of optimal solutions to examples E1-E4, with different weak bond strengths

# created Jul 10 2025

# SP iterator -- see https://stackoverflow.com/questions/57445750/iteration-index-in-a-julia-for-loop

# ****REMEMBER****: must input NEGATIVE delta values. 

# TO DO:
# - 


using Plots,  LaTeXStrings, DifferentialEquations

include("E0.jl")
include("E1.jl")
include("E2.jl")
include("E3.jl")
include("E4.jl")
include("data_weakvary.jl")


# set up plots
default(fontfamily = "Computer Modern")
default(linewidth=2,titlefontsize=16,legendfontsize=12,legendtitlefontsize=14,
        guidefontsize=14, tickfontsize=12)

my_dpi = 300

# set up subplots
l = @layout [grid(1,4) a{0.15w}]
plt = plot(
    layout=l,
    link=:y,
    size=(1300,250), #(round(Int, 17.8 / 2.54 * my_dpi), round(Int, 3 / 2.54 * my_dpi)),
    #dpi=my_dpi,
    margin=0Plots.mm,
    bottom_margin=12Plots.mm,
    left_margin=12Plots.mm
)   
lstrings = ["10^3"; "10^4"; "10^5"]


# add title/axis labels
for i in 1:4
    plot!(
        xlabel=L"\delta_w",
        ylabel=(i == 1 ? "yield" : ""),
        yticks=(i == 1 ? :auto : false),
        ylimits=(0, 1),
        title="E" * string(i),
        subplot=i
    )
end


# get colour palettes
copt1 = palette(:PuBu_4,rev=false)
copt2 = palette(:OrRd_4,rev=false)
ceq = palette(:Greens_4)


# set up legend
plot!(plt, showaxis = false, grid = false, legend=(-0.2,0.7),  subplot=5) #legend=(-0.3,0.7)
for it=1:3
    plot!(plt,[],label=latexstring(lstrings[it])*" ST",color=copt1[it+1], subplot=5)
end
for it=1:3
    plot!(plt,[],label=latexstring(lstrings[it])*" EQ",color=ceq[it+1],alpha = 0.8,subplot=5)
end
#plot!(plt,leg_title = "T", legend=:inside,subplot=5)


# parameters
tmaxlist = [1000.0,10000.0,100000.0]    # max time
ds = 15.0   # strong bond
dwlist = collect(0:0.5:ds)  # weak bonds
dens = 0.05;   # density



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
        ymax = 0.8656837320424821;
        f! = E1!
    elseif whichcase==2
        nv = 7
        nblocks = 95
        del1=[0,1,1,1,1,1,0,1,0,0,0,0]
        del2=[1,0,0,0,0,0,1,0,1,1,1,1]
        ymax = 0.9194026669911366
        f! = E2!
    elseif whichcase==3
        nv = 8
        nblocks = 167
        del1=[0,1,0,1,0,1,1,0,0,1,1,1]
        del2=[1,0,1,0,1,0,0,1,1,0,0,0]
        ymax = 0.8938384955289904
        f! = E3!
    elseif whichcase ==4
        nv = 6
        nblocks = 63
        del1=[1,0,1,1,0,0,0,0,0,1,0,0,0,1,0]
        del2=[0,1,0,0,1,1,1,1,1,0,1,1,1,0,1]
        ymax = 0.9416962259076676
        f! = E4!
    end

    # initial concentration
    c0 = zeros(nblocks);
    c0[1:nv] .= dens/nv;

    # set up ode problem
    prob = ODEProblem(f!,c0,(0,10),-del1)

    # solve it
    odesolver = AutoTsit5(Rosenbrock23())  # this one works way better than Tsit5()
    #@time sol = solve(prob,odesolver)  
    
    #plot(sol.t,sol[nblocks,:]/dens*nv, linewidth = 2, xaxis="t", yaxis = "yield") 
    #plot(sol)

    # define function to allow varying parameters
    function sol_params(deltas,t)
        newproblem = remake(prob,p = -deltas, tspan = (0.0,t))
        sol = solve(newproblem, odesolver, save_on = false)
        yield = sol[end,end] #Final concentration of target structure
        return yield/(dens/nv)
    end


    # plot max yield previously obtained 
    plot!(plt,[0,ds],[ymax,ymax],linestyle=:dash,linewidth=1,linecolor=:black,label="",subplot=sp)
    #plot!(plt,[0,ds],[ymax,ymax],linestyle=:dash,linewidth=1,linecolor=:black,label="Optimal "*L"(10^4)",subplot=sp)
    
  
    # loop through and load data and plot it
    # (option to generate it too, if uncomment some code)
    for (it,tmax) in enumerate(tmaxlist)

        println("tmax = ",tmax)

        # uncomment this, to generate data
        # # get equal energies
        # println("equal bonds")
        # yieldEq = zeros(size(dwlist))
        # for i in 1:length(dwlist)
        #     dw = dwlist[i]
        #     deltas = dw*del1 + dw*del2   # bond strengths
        #     yieldEq[i] = sol_params(deltas,tmax)
        # end
        # println(yieldEq) # print data to output, so we can save it for later

        # load data that we generated previously
        #yieldEq = ydata(whichcase, tmax,1)

        # plot equal energies
        #plot!(plt,dwlist,yieldEq,label="",linestyle=:dash,subplot=sp)


        # uncomment this, to generate data
        # # solve for different weak bonds
        # println("strong/weak bonds")
        # yield = zeros(size(dwlist))
        # for i in 1:length(dwlist)
        #     dw = dwlist[i]
        #     deltas = ds*del1 + dw*del2   # bond strengths
        #     yield[i] = sol_params(deltas,tmax)
        # end
        # println(yield) # print data to output, so we can save it for later


        # load data that we generated previously
        yield = ydata(whichcase, tmax, 0)

        # plot
        plot!(plt,dwlist,yield,label=latexstring(lstrings[it]),legend=:topright,color=copt1[it+1],subplot=sp)
        #plot!(plt,dwlist,yield,label=string(Int64(tmax)),
        #    legend=:topright,title="E"*string(whichcase),subplot=sp)

        # println("whichcase = ",whichcase)
        # for j=1:length(yield)
        #     println(dwlist[j],", ",yield[j])
        # end

    end  # loop through tmax

    # plot equal energies
    for (it,tmax) in enumerate(tmaxlist)
        yieldEq = ydata(whichcase, tmax,1)
        plot!(plt,dwlist,yieldEq,label=L"10^4"*"(equal)",color=ceq[it+1],alpha = 0.8,subplot=sp)
    end

    # # plot equal energies for tmax=10000
    # yieldEq = ydata(whichcase, 10000,1)
    # plot!(plt,dwlist,yieldEq,label=L"10^4"*"(equal)",linestyle=:dash,color=cdefault[4],alpha = 0.6,subplot=sp)
    
    #if whichcase != 4
        plot!(plt,legend=false,subplot=sp)
    #end

end  # loop through whichcase

display(plt)

# save (adjust size until font sizes are large enough when saved)
savefig("./Figure3.pdf")