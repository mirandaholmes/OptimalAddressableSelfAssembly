#= 
    plot output from compute_dwopt.jl

    created nov 1 2025


=#

using DelimitedFiles, Plots, LaTeXStrings, ColorSchemes, Polynomials

ifsave = false

times = 10.0.^collect(3:0.25:5); #10.0.^collect(3:0.25:5);


# set up plot
default(fontfamily = "Computer Modern")
default(linewidth=2,titlefontsize=16,legendfontsize=10,legendtitlefontsize =12,guidefontsize = 14, tickfontsize=12)
plsize = (600,450)
plt = plot(layout=(2,2),size=plsize,xscale=:log10,xlabel="T",ylabel=L"\delta_w",)#,bottom_margin=8Plots.mm,left_margin=5Plots.mm)


for wc in [1,2,3,4]

    datafile = "data_dwE"*string(wc)*".txt"
    data = readdlm(datafile)

    ds = data[:,1]
    dw = data[:,2:end]

    
    plot!(plt[wc],times,dw',label=["15" "20" "25"],
            title="E"*string(wc),palette=:Dark2_3,marker=:circle,legendtitle=L"\delta_s")

    # fit a line through ds=15 data
    istart = 2;
    xs = log.(times[istart:end]);
    ys = dw[1,istart:end];
    p = fit(xs,ys,1);
    println("wc=",wc,", slope = ",p[1])
    

end

display(plt)

if ifsave
    savefig(plt,"dwopt.pdf")
end
