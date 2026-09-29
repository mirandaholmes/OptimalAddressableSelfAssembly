#=
    plot yields for all the spanning trees, at given values of delta

    created Oct 29 2025


    format of saved data in datafiles:
        [ST#  ds  dw  yields(times)]
        ... (1 row per ST) ..

    times = 10.0.^collect(3:0.25:5) = [1000.0  1778.28  3162.28  5623.41  10000.0  17782.8  31622.8  56234.1  100000.0]'
=#

using DelimitedFiles, Plots, LaTeXStrings, ColorSchemes

# parameters
whichcase = [1,2,3,4]
ds = -15.0
dw = -0.0
tind = 5     # which time index to plot (5=10^4)
times = 10.0.^collect(3:0.25:5)

# set up plot
plt = plot(ylims=(0,1))
cd = palette(:default)

for wc in whichcase
    # read in data
    datafile = "dataYieldsE"*string(wc)*".txt"   # save data here
    data = readdlm(datafile)

    # extract data corresponding to desired parameters
    rowinds = findall((data[:,2] .== ds) .& (data[:,3] .== dw))
    yields = data[rowinds,4:end];
    ypl = yields[:,tind]
    sort!(ypl,rev=true)

    # plot
    #scatter!(ypl,label="E"*string(wc))

    # plot all trajectories instead
    plot!(times, yields[:,:]',xscale=:log10,color=cd[wc], alpha = 0.05, linewidth=2, label=:none)
end


display(plt)

