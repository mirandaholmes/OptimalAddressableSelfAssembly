#=
    plot figures for paper showing correlations
        (a) yield vs bad volumes, all examples
        (b) bad volumes vs S, all examples
        (c) S vs W, all examples
        (d) yield vs S, all examples

    Note: here S=rho-index

    created Oct 30 2025


    format of saved data in datafiles:
        [ST#  ds  dw  yields(times)]
        ... (1 row per ST) ..

    times = 10.0.^collect(3:0.25:5) = [1000.0  1778.28  3162.28  5623.41  10000.0  17782.8  31622.8  56234.1  100000.0]'

wc = 1
 corr(badvol, yield) = -0.973345326792357
 corr(badvol, S) = -0.9051338204208625
 corr(S,W) = -0.930787635327742
wc = 2
 corr(badvol, yield) = -0.9770356646290673
 corr(badvol, S) = -0.748789071623088
 corr(S,W) = -0.9801238966718406
wc = 3
 corr(badvol, yield) = -0.9868723070098228
 corr(badvol, S) = -0.8748641588736513
 corr(S,W) = -0.9175292645280151
wc = 4
 corr(badvol, yield) = -0.9483847270386749
 corr(badvol, S) = -0.9161528051630061
 corr(S,W) = -0.9921835635830106
=#

using DelimitedFiles, Plots, LaTeXStrings, ColorSchemes, Statistics

include("HelpersInvariants.jl")
include("HelpersSubTrees.jl")

# parameters
ifsave = false
ds = -15.0
dwlist = -[5.0,3.5,3.5,2.5]; 
tinds = [1,5,9]     # which time index to plot yield in plot#3 (5=10^4)
times = 10.0.^collect(3:0.25:5)
nt = length(times)

# savefiles
savefile1 = "corr_yieldBadvol.pdf"
savefile2 = "corr_badvolRho.pdf"
savefile3 = "corrWrho.pdf"
savefile4 = "corrYieldRho.pdf"


# set up plots -- 3 of them
default(fontfamily = "Computer Modern")
default(linewidth=2,titlefontsize=16,legendfontsize=10,legendtitlefontsize =12,guidefontsize = 14, tickfontsize=12)
plsize = (800,600)
plt1 = plot(layout=(2,2),size=plsize,xlabel="Off-pathway volume (%)",ylabel="Yield",ylims=(0,1.05))#,bottom_margin=8Plots.mm,left_margin=5Plots.mm)
plt2 = plot(layout=(2,2),size=plsize,xlabel="Off-pathway volume (%)",ylabel=L"\rho")#,bottom_margin=8Plots.mm,left_margin=5Plots.mm)
plt3 = plot(layout=(2,2),size=plsize,xlabel=L"\rho",ylabel="W")#,bottom_margin=8Plots.mm,left_margin=5Plots.mm)
plt4 = plot(layout=(2,2),size=plsize,xlabel=L"\rho",ylabel="Yield")#,bottom_margin=8Plots.mm,left_margin=5Plots.mm)


# loop through examples and load data
for wc in [1,2,3,4]

    dw = dwlist[wc]

    # compute graph invariants
    deltafile = "dataTreesE"*string(wc)*".txt"
    deltas = readdlm(deltafile,skipstart=3)
    numtrees = size(deltas)[1]
    S,W = getInvariants(wc);

    # read in yield data
    datafile = "dataYieldsE"*string(wc)*".txt"   # save data here
    data = readdlm(datafile)

    # extract data corresponding to desired parameters
    rowinds = findall((data[:,2] .== ds) .& (data[:,3] .== dw))
    yields = data[rowinds,4:4+nt-1];
    badvols = data[rowinds,4+nt:end]
    ypl = yields[:,tinds]
    badpl = badvols[:,tinds]


    # ------------------------------------
    #   Plot #1 -- yield vs bad volumes
    # ------------------------------------
    scatter!(plt1[wc],badpl*100,ypl, color_palette =:Dark2_3, markershape=[:diamond :square :circle],
            title="E"*string(wc)*", "*L"\delta_w="*string(-dw),label=[L"10^3" L"10^4" L"10^5"],legendtitle="T")

    println("wc = ",wc)
    println(" corr(badvol, yield) = ",cor(badpl[:,2],ypl[:,2]))

    # ------------------------------------
    #   Plot #2 -- S vs bad volumes
    # ------------------------------------
    scatter!(plt2[wc],badpl[:,2]*100,S, c=palette(:Dark2_3)[2], markershape=[ :square :circle],
            title="E"*string(wc)*", "*L"\delta_w="*string(-dw),label=[L"10^4" L"10^5"],legendtitle="T")

    println(" corr(badvol, S) = ",cor(badpl[:,2],S))

    # ------------------------------------
    #   Plot #3 -- W vs S
    # ------------------------------------
    scatter!(plt3[wc],S,W, c=palette(:Dark2_4)[4], markershape=:circle,
            title="E"*string(wc),label=:none)

    println(" corr(S,W) = ",cor(S,W))

    # ------------------------------------
    #   Plot #4 -- yield vs S
    # ------------------------------------
    scatter!(plt4[wc],S,ypl, color_palette =:Dark2_3, markershape=[:diamond :square :circle],
            title="E"*string(wc)*", "*L"\delta_w="*string(-dw),label=[L"10^3" L"10^4" L"10^5"],legendtitle="T")

    println(" corr(S,Yield) = ",cor(S,ypl[:,2]))
    println(" corr(W,Yield) = ",cor(W,ypl[:,2]))
end


display(plt4)

if ifsave
    savefig(plt1,savefile1)
    savefig(plt2,savefile2)
    savefig(plt3,savefile3)
    savefig(plt4,savefile4)
end
