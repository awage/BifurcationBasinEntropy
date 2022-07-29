using JLD2
using CairoMakie
using LaTeXStrings
using DynamicalSystems
@load "results_henon_res_500.jld2"
a = range(0., 3.98, length = 1000)
ν = 0.01
@load "diag_bif.jld2"
fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = "y", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20)
ax2 = Axis(fig[2,1], ylabel = "Sb", xlabel = "F", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)

#ax = Axis(fig[1,1], ylabel = "y", yticklabelsvisible = false, xticklabelsvisible = false, bottomspinevisible = false, topspinevisible = false, leftspinevisible = false, rightspinevisible = false, xgridvisible = false, ygridvisible = false)
#ax2 = Axis(fig[2,1], ylabel = "Sb", xlabel = "F")

ind = findall(P[:,1] .≤ 4)
scatter!(ax, P[ind,1],P[ind,2], markersize = 0.7, color = :black, rasterize = 4)

ind = findall(a .≤ 4)

lines!(ax2, a[ind], vec(Sb[ind]))
#lp.rasterize = 4

save("diag_bif_henon.svg",fig)
