#=
    plot figures for paper: 
        (a) yields vs time for dw=0  (all trees)
        (b) yields vs time for dw > 0  (all trees)
        (c) yield vs N for dw > 0, T=10^4

    created Oct 29 2025


    format of saved data in datafiles:
        [ST#  ds  dw  yields(times)]
        ... (1 row per ST) ..

    times = 10.0.^collect(3:0.25:5) = [1000.0  1778.28  3162.28  5623.41  10000.0  17782.8  31622.8  56234.1  100000.0]'

=#

using DelimitedFiles, Plots, LaTeXStrings, ColorSchemes, Statistics

include("HelpersInvariants.jl")
include("HelpersSubTrees.jl")

# parameters
wc = 4
ifsave = false
savefile = "yieldsE"*string(wc)*".pdf"
ds = -15.0
dwlist = -[5.0,3.5,3.5,2.5]; dw = dwlist[wc]
    # - e1: 5,5.5
    # - e2: 3,3.5,4
    # - e3: 2.5,3,3.5
    # - e4: 2 = big decrease, could also try 1.5
tind = [1,5,9]     # which time index to plot yield in plot#3 (5=10^4)
times = 10.0.^collect(3:0.25:5)
nt = length(times)



# compute graph invariants
deltafile = "dataTreesE"*string(wc)*".txt"
deltas = readdlm(deltafile,skipstart=3)
numtrees = size(deltas)[1]
S,W = getInvariants(wc);
indw = sortperm(W);   # order from smallest to largest


# read in yield data
datafile = "dataYieldsE"*string(wc)*".txt"   # save data here
data = readdlm(datafile)


# set up plot
default(fontfamily = "Computer Modern")
default(linewidth=2,titlefontsize=12,legendfontsize=8,legendtitlefontsize =10,guidefontsize = 14, tickfontsize=12)
#l = @layout [(1,3) a{0.1w}]
plt = plot(layout=(1,3),size=(800, 200),bottom_margin=8Plots.mm,left_margin=5Plots.mm)

# colour palette
cdef = palette(:default)   #[cdef[1] ; cdef[1]]
wc < 4 ? cd = cgrad([cdef[1] ; cdef[1]], numtrees, categorical = true) : # colour scheme for E4
         cd = palette(:Dark2_6)    # colour scheme for E1,E2,E3
 
# alpha values of lines and markers, for each case
aline = [0.1 0.1 0.1 1]
amarker = [1 1 1 1]



# extract data corresponding to desired parameters
rowinds0 = findall((data[:,2] .== ds) .& (data[:,3] .== 0))
rowinds1 = findall((data[:,2] .== ds) .& (data[:,3] .== dw))


# ------------------------------------
#   Plot #1 -- yield vs time, dw=0
# ------------------------------------

ypl0 = data[rowinds0,4:4+nt-1]'
if wc == 4
    plot!(plt,times, ypl0[:,indw], color_palette=cd, alpha=aline[wc],  label=["Tree 1" "Tree 2" "Tree 3" "Tree 4" "Tree 5" "Tree 6"], 
           subplot=1)
else
    plot!(plt,times, ypl0[:,indw], color_palette=cd, alpha=aline[wc],  label=:none, subplot=1)
end
plot!(plt,xlabel=L"T", ylabel="Yield", xscale=:log10, ylims=(0,1.05), title="E"*string(wc)*", "*L"\delta_w=0",subplot=1)


# ------------------------------------
#   Plot #2 -- yield vs time, dw>0
# ------------------------------------

ypl1 = data[rowinds1,4:4+nt-1]'
plot!(plt,times, ypl1[:,indw], xscale=:log10,color_palette=cd, alpha=aline[wc],label=:none, subplot=2)
plot!(plt,xlabel=L"T", ylims=(0,1.05), title=L"\delta_w="*latexstring(string(-dw)),subplot=2)
plot!(plt, ylabel="", yformatter=_->"", left_margin=0Plots.mm, subplot=2)

# ------------------------------------
#   Plot #3 -- yield vs N(T), dw>0
# ------------------------------------
yields = ypl1';
ypl = yields[:,tind]
if wc < 4
    scatter!(plt,W,reverse(ypl,dims=2),alpha=amarker[wc], label=:none, seriescolor=[cdef[3] cdef[2] cdef[1]], subplot=3)
    scatter!(plt,[W[4]],2*[1 1 1],seriescolor=[cdef[1] cdef[2] cdef[3]],label = [L"10^3" L"10^4" L"10^5"], subplot=3)
else
    # plot markers, with deisgnated colours
    for j=1:6
        scatter!(plt[3],[W[j]],[ypl[j,3]],label=:none,c=cd[j],marker=:circle)
        scatter!(plt[3],[W[j]],[ypl[j,2]],label=:none,c=cd[j],marker=:square)
        scatter!(plt[3],[W[j]],[ypl[j,1]],label=:none,c=cd[j],marker=:diamond)
    end
    # make legend
    scatter!(plt[3],[W[4]],2*[1],c=:grey,label = L"10^3", marker=:diamond)
    scatter!(plt[3],[W[4]],2*[1],c=:grey,label = L"10^4", marker=:square)
    scatter!(plt[3],[W[4]],2*[1],c=:grey,label = L"10^5", marker=:circle)
end
plot!(plt,xlabel=L"\mathcal{W}",  ylims=(0,1.05), title=L"\delta_w="*latexstring(string(-dw)),
         legend=:outerright, legendtitle=L"T", subplot=3)
plot!(plt,ylabel="", yformatter=_->"", left_margin=0Plots.mm, subplot=3)


# compute correlations
println("T1")
println("wc = ",wc,", corr(W,yield) = ",cor(W,ypl[:,1]))
println("wc = ",wc,", corr(S,yield) = ",cor(S,ypl[:,1]))
println("T2")
println("wc = ",wc,", corr(W,yield) = ",cor(W,ypl[:,2]))
println("wc = ",wc,", corr(S,yield) = ",cor(S,ypl[:,2]))
println("T3")
println("wc = ",wc,", corr(W,yield) = ",cor(W,ypl[:,3]))
println("wc = ",wc,", corr(S,yield) = ",cor(S,ypl[:,3]))

display(plt)

if ifsave
    savefig(savefile)
end
