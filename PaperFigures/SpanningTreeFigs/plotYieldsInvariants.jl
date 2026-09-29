#=
    plot yields for all the spanning trees vs invariants

    created Oct 29 2025


    format of saved data in datafiles:
        [ST#  ds  dw  yields(times)]
        ... (1 row per ST) ..

    times = 10.0.^collect(3:0.25:5) = [1000.0  1778.28  3162.28  5623.41  10000.0  17782.8  31622.8  56234.1  100000.0]'
=#

using DelimitedFiles, Plots, LaTeXStrings, ColorSchemes

include("HelpersInvariants.jl")
include("HelpersSubTrees.jl")

# parameters
whichcase = 1#[1,2,3,4]
ds = -15.0
dw = -3.0
tinds = [1,5,9]     # which time index to plot (5=10^4)
times = 10.0.^collect(3:0.25:5)
nt = length(times);

# set up plot
plt = plot()
#cd = palette(:default)


for wc in whichcase
    #wc=4

    # compute graph invariant
    N,W = getInvariants(wc);

    # read in yield data
    datafile = "dataYieldsE"*string(wc)*".txt"   # save data here
    data = readdlm(datafile)

    # extract data corresponding to desired parameters
    rowinds = findall((data[:,2] .== ds) .& (data[:,3] .== dw))
    yields = data[rowinds,4:4+nt-1];
    badvols = data[rowinds,4+nt:end]
    ypl = yields[:,tinds]
    badpl = badvols[:,tinds]
    #ypl = yields[:,tind]

    # plot
    #scatter!(N,ypl,alpha=1,label="E"*string(wc),ylims=(0,1),xlabel="N",ylabel="yield")
    scatter!(W,ypl,alpha=1,label="E"*string(wc),ylims=(0,1),xlabel="W",ylabel="yield")
    #scatter!(badpl*100,ypl,xlabel="vol",ylabel="Yield",ylims=(0,1))
    #scatter!(N,badpl*100,xlabel="N",ylabel="Vol")
    #scatter!(N,W,xlabel="N",ylabel="W")
end
display(plt)

