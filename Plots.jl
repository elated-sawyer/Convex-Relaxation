#Performance profiles for Mosek, Cplex and Gurobi.
using DelimitedFiles
using BenchmarkProfiles, Plots
dat = readdlm( "Results for Solvers.in" )
for i in 1:1:99
    for j in 1:1:3
        if dat[i,j]>=360
            dat[i,j] = Inf
        end
    end
end
p = performance_profile(PlotsBackend(), dat, 
    ["Mosek", "Cplex", "Gurobi"], 
    legend = :topright,
    ylims =(0,0.35),
    #xaxis = 
    #scale = :none,
    )
savefig(p,"ComSol1.pdf")

#DNN-SOCP detail process
using CSV
using DataFrames
using LazyStack
using Plots
df = CSV.read("Records for Iterations.csv", DataFrame)

#Iteration process of SOCP-based cutting plane method.
aa = [4,19,43]
dim = ["Dimension 30","Dimension 40","Dimension 50"]
plot_array = Any[]
for i in [1,2,3]
    k = aa[i]
    k = (k-1)*7+1
    k_n = k
    z_name = [SubString(df[1,k_n],5),SubString(df[1,k_n+7],5),SubString(df[1,k_n+14],5),SubString(df[1,k_n+21],5),SubString(df[1,k_n+28],5)]
    k_x = k+1
    x_time = rstack(df[:,k_x], df[:,k_x+7], df[:,k_x+14], df[:,k_x+21], df[:,k_x+28]; fill=missing)
    k_y = k+3
    y_gap = rstack(df[:,k_y], df[:,k_y+7], df[:,k_y+14], df[:,k_y+21], df[:,k_y+28]; fill=missing)
    p = plot(x_time, y_gap,
            #markershape = :hexagon,
            markeralpha = 0.6,
            markersize = 5,
            markerstrokewidth = 0.1,
            markerstrokealpha = 0.1,
            #markerstrokecolor = :black,
            #markerstrokestyle = :dot,
            seriestype = :scatter,
            label = permutedims(z_name),
            legendfontsize = 5,
            legend = :topright,
            title = dim[i],
            titlelocation = :left,
            titlefontsize = 9,
            legendfont = 6,
            yticks=-0.2:0.2:1, ylims =(-0.2,1))
    push!(plot_array,p)
end
#for dimension 60
k = 52
k = (k-1)*7+1
k_n = k
z_name = [SubString(df[1,k_n],5),SubString(df[1,k_n+7],5),SubString(df[1,k_n+14],5)]
k_x = k+4
x_time = rstack(df[:,k_x], df[:,k_x+7], df[:,k_x+14]; fill=missing)
k_y = k+3
y_gap = rstack(df[:,k_y], df[:,k_y+7], df[:,k_y+14]; fill=missing)
p = plot(x_time, y_gap,
        #markershape = :hexagon,
        markeralpha = 0.6,
        markersize = 5,
        markerstrokewidth = 0.1,
        markerstrokealpha = 0.1,
        #markerstrokecolor = :black,
        #markerstrokestyle = :dot,
        seriestype = :scatter,
        label = permutedims(z_name),
        legendfontsize = 5,
        legend = :topright,
        title = "Dimension 60",
        titlelocation = :left,
        titlefontsize = 9,
        legendfont = 6,
        yticks=-0.2:0.2:1, ylims =(-0.2,1))
push!(plot_array,p)
plot(plot_array..., 
    xlabel="Time(s)", ylabel="Relative gap (%)",
    #xtickfontsize=18,
    xguidefontsize=9,yguidefontsize=9,
    layout=(2,2))
#plot!(size=(700,400))
savefig("Low1.pdf")

#Iteration process of SOCP-based cutting plane method for high dimensional instances.
aa = [64,73,82,91]
dim = ["Dimension 80","Dimension 90","Dimension 100","Dimension 125"]
plot_array = Any[]
for i in [1,2,3,4]
    k = aa[i]
    k = (k-1)*7+1
    k_n = k
    z_name = [SubString(df[1,k_n],5),SubString(df[1,k_n+7],5),SubString(df[1,k_n+14],5),SubString(df[1,k_n+21],5),SubString(df[1,k_n+28],5)]
    k_x = k+4
    x_time = rstack(df[:,k_x], df[:,k_x+7], df[:,k_x+14], df[:,k_x+21], df[:,k_x+28]; fill=missing)
    k_y = k+3
    y_gap = rstack(df[:,k_y], df[:,k_y+7], df[:,k_y+14], df[:,k_y+21], df[:,k_y+28]; fill=missing)
    p = plot(x_time, y_gap,
            #markershape = :hexagon,
            markeralpha = 0.6,
            markersize = 5,
            markerstrokewidth = 0.1,
            markerstrokealpha = 0.1,
            #markerstrokecolor = :black,
            #markerstrokestyle = :dot,
            seriestype = :scatter,
            label = permutedims(z_name),
            legendfontsize = 5,
            legend = :topright,
            title = dim[i],
            titlelocation = :left,
            titlefontsize = 9,
            legendfont = 6,
            #ylims =(-0.2,1),
            #yticks=-0.2:0.2:1
            )
    push!(plot_array,p)
end
plot(plot_array..., 
    xlabel="Time(s)", ylabel="Relative gap (%)",
    #xtickfontsize=18,
    xguidefontsize=9,yguidefontsize=9,
    layout=(2,2))
#plot!(size=(700,400))
savefig("HigSca1.pdf")

#The iterative process of SOCP-based cutting plane method for high dimensional instance.
aa = [64,73,82,91]
dim = ["Dimension 80","Dimension 90","Dimension 100","Dimension 125"]
bb = ["Instance","Computing time for each iteration","Number of added constraints","Relative gap to optimal value (%)","Cumulative computing time","The number of iterations","Total number of constraints"]
plot_array = Any[]
for m in [3,4,2,5]
    k = 64
    k = (k-1)*7+1
    k_n = k
    z_name = [SubString(df[1,k_n],5),SubString(df[1,k_n+7],5),SubString(df[1,k_n+14],5),SubString(df[1,k_n+21],5),SubString(df[1,k_n+28],5)]
    k_x = k+5
    x_time = rstack(df[:,k_x], df[:,k_x+7], df[:,k_x+14], df[:,k_x+21], df[:,k_x+28]; fill=missing)
    k_y = k+m-1
    y_gap = rstack(df[:,k_y], df[:,k_y+7], df[:,k_y+14], df[:,k_y+21], df[:,k_y+28]; fill=missing)
    if m == 4
        p = plot(x_time, y_gap, #w_tim,
        #markershape = :hexagon,
        markeralpha = 0.6,
        markersize = 5,
        markerstrokewidth = 0.1,
        markerstrokealpha = 0.1,
        #markerstrokecolor = :black,
        #markerstrokestyle = :dot,
        #seriestype = :scatter,
        seriestype = :steppre,
        label = permutedims(z_name),
        legendfontsize = 5,
        legend = :topright,
        title = bb[m],
        titlelocation = :left,
        titlefontsize = 9,
        legendfont = 6,
        #ylims =(-0.2,1),
        #yticks=-0.2:0.2:1
        )
        push!(plot_array,p)
    elseif m == 5
        p = plot(x_time, y_gap, #w_tim,
        #markershape = :hexagon,
        markeralpha = 0.6,
        markersize = 5,
        markerstrokewidth = 0.1,
        markerstrokealpha = 0.1,
        #markerstrokecolor = :black,
        #markerstrokestyle = :dot,
        #seriestype = :scatter,
        seriestype = :steppre,
        label = permutedims(z_name),
        legendfontsize = 5,
        legend = :topleft,
        title = bb[m],
        titlelocation = :left,
        titlefontsize = 9,
        legendfont = 6,
        #ylims =(-0.2,1),
        #yticks=-0.2:0.2:1
        )
        push!(plot_array,p) 
    else
        p = plot(x_time, y_gap, #w_tim,
        #markershape = :hexagon,
        markeralpha = 0.6,
        markersize = 5,
        markerstrokewidth = 0.1,
        markerstrokealpha = 0.1,
        #markerstrokecolor = :black,
        #markerstrokestyle = :dot,
        #seriestype = :scatter,
        seriestype = :steppre,
        label = permutedims(z_name),
        legendfontsize = 5,
        legend = :bottomright,
        title = bb[m],
        titlelocation = :left,
        titlefontsize = 9,
        legendfont = 6,
        #ylims =(-0.2,1),
        #yticks=-0.2:0.2:1
        )
        push!(plot_array,p)
    end
end
plot(plot_array..., 
    xlabel="The number of iterations", #ylabel="Gap(%)",
    #xtickfontsize=18,
    xguidefontsize=9,yguidefontsize=9,
    layout=(2,2))
#plot!(size=(700,400))
savefig("AppPro2.pdf")

# The change of computing time for single iteration under different dimensions.
aa = [55,64,73,82]
dim = ["Dimension 70","Dimension 80","Dimension 90","Dimension 100"]
plot_array = Any[]
for i in [1,2,3,4]
    k = aa[i]
    k = (k-1)*7+1
    k_n = k
    z_name = [SubString(df[1,k_n],5),SubString(df[1,k_n+7],5),SubString(df[1,k_n+14],5),SubString(df[1,k_n+21],5),SubString(df[1,k_n+28],5)]
    k_x = k+5
    x_time = rstack(df[:,k_x], df[:,k_x+7], df[:,k_x+14], df[:,k_x+21], df[:,k_x+28]; fill=missing)
    k_y = k+1
    y_gap = rstack(df[:,k_y], df[:,k_y+7], df[:,k_y+14], df[:,k_y+21], df[:,k_y+28]; fill=missing)
    p = plot(x_time, y_gap,
            #markershape = :hexagon,
            markeralpha = 0.6,
            markersize = 5,
            markerstrokewidth = 0.1,
            markerstrokealpha = 0.1,
            #markerstrokecolor = :black,
            #markerstrokestyle = :dot,
            seriestype = :scatter,
            label = permutedims(z_name),
            legendfontsize = 5,
            legend = :topleft,
            title = dim[i],
            titlelocation = :left,
            titlefontsize = 9,
            legendfont = 6,
            #ylims =(-0.2,1),
            #yticks=-0.2:0.2:1
            )
    push!(plot_array,p)
end
plot(plot_array..., 
    xlabel="The number of iterations", ylabel="Time(s)",
    #xtickfontsize=18,
    xguidefontsize=9,yguidefontsize=9,
    layout=(2,2))
#plot!(size=(700,400))
savefig("AppPro3.pdf")



