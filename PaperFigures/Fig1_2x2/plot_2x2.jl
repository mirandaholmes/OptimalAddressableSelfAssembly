#=
Plot Fig. 1, studying optimization for a 2x2 square

=======   
 TO DO:
=======
- 
=#

using Plots, DelimitedFiles, LaTeXStrings, Colors, Polynomials
include("eqyield.jl")


ifsave = 0  # set to 1 to save data

# data file and save file
datafile = joinpath(@__DIR__,"data.txt")
savefile = joinpath(@__DIR__,"Fig1_2x2.pdf")

# read data
# format: density, time, c4, gradient, optimal deltas
data = readdlm(datafile,',')
times = data[:,2]
yields = data[:,3]
dels = data[:,5:end];

# sort delta-values
sort!(dels,dims=2)

# colours for plotting
colorred = RGB(239/255,138/255,98/255)
colorblue = RGB(103/255,169/255,207/255)
colorgreen = RGB(141/255,211/255,199/255)
colorpurple = RGB(190/255,186/255,218/255)

# set up plot
plt=plot(xlabel=L"\textrm{Time} \;\; T",ylabel="Interaction Energy",legend=:bottomright)
plot!(guidefont=(16,"Computer Modern"),xtickfont=font(15,"Computer Modern"),ytickfont=font(15,"Computer Modern"),
legendfont=font(15,"Computer Modern"),titlefont=font(17,"Computer Modern"),
ylim=(-0.5,16))



# Plot deltas
plot!(times,dels[:,1],xaxis=:log10,linewidth=3,markershape=:xcross,markersize=7,label=L"\delta_{34}",color=colorred)
plot!(times,dels[:,2],xaxis=:log10,linewidth=3,markershape=:utriangle,markersize=9,label=L"\delta_{12}",color=colorpurple)
plot!(times,dels[:,3],xaxis=:log10,linewidth=3,markershape=:rect,markersize=8,label=L"\delta_{13}",color=colorblue)
plot!(times,dels[:,4],xaxis=:log10,linewidth=3,markershape=:circle,markersize=7,label=L"\delta_{24}",color=colorgreen)

# fit a line through delta_34
istart = 7;
xs = log.(times[istart:end]);
ys = dels[istart:end,1];
p = fit(xs,ys,1);
println("slope = ",p[1])


# plot yield, on a new axis
axy = twinx(plt)
plot!(axy, times,yields, color=:grey,linewidth=3,alpha=0.6,
      xaxis=:log10, label="Yield",ylabel="Yield",ylim=(0,1.0),
      legend=:none,
      ytickfont=font(15,"Computer Modern",:grey),
      legendfont=font(15,"Computer Modern"),
      guidefont=(16,"Computer Modern",:grey))



# calculate equilibrium yield, for given delta-values
eqyields = similar(yields)
for it=1:length(times)
    dw=dels[it,1];
    ds=dels[it,2];
    eqyields[it]=ceq(ds,dw,0.05/4);
    #println("dw=",dw,",ds=",ds,",eq=",eqyields[it])
end
plot!(axy,times,eqyields,color=:grey,linestyle=:dash,linewidth=2,alpha=0.6,
      legend=:none)



# display plot
display(plt)

if ifsave==1
    savefig(plt,savefile)
end
