using JLD2
using CairoMakie
using LaTeXStrings
using DynamicalSystems
@load "diag_bif.jld2"
fig = Figure(resolution = (800, 600))
ax2 = Axis(fig[1,1], ylabel = "Sb", xlabel = "a", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
ind = findall(0.97 .≤ a .≤ 1.02)
lines!(ax2, a[ind], vec(Sb[ind]))
save("zoom_bif_1.svg",fig)

fig = Figure(resolution = (800, 600))
ax2 = Axis(fig[1,1], ylabel = "Sb", xlabel = "a", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
ind = findall(2.33 .≤ a .≤ 2.50)
lines!(ax2, a[ind], vec(Sb[ind]))
save("zoom_bif_2.svg",fig)


fig = Figure(resolution = (800, 1000))
ax2 = Axis(fig[1,1], ylabel = "Sb", xlabel = "a", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
ind = findall(2.5 .≤ a .≤ 2.80)
lines!(ax2, a[ind], vec(Sb[ind]))
ax3 = Axis(fig[2,1], ylabel = "Na", xlabel = "a", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
Na = [length(keys(a)) for a in att]
lines!(ax3, a[ind], vec(Na[ind]))
save("zoom_Na_1.svg",fig)



