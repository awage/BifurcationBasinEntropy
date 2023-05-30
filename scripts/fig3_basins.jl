# Script to produce the Figure 3 of the paper 
# Using basin entropy to explore bifurcations
# MIT License
# Copyright (c) 2023 Alexandre Wagemkakers
# You should have received a copy of the license with this software. 
using DrWatson
@quickactivate
using Attractors
using CairoMakie
using LaTeXStrings

function pitchfork!(dz, z, p, t)
    x, y = z
    dz[1] = p*x - x^3 
    dz[2] = 0.1*y 
    return
end

function _get_pitchfork(μ, xg, yg)
    ds = DeterministicIteratedMap(pitchfork!, [1.0, 0.0], μ)
    xgg = range(-2, 2, length = 5000)
    ygg = range(-1, 1, length = 5000)
    grid = (xgg, ygg)
    mapper = AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 10000,
            mx_chk_loc_att = 10000,
            mx_chk_att = 2)
    grid = (xg, yg)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end

function compute_Sb_fig(μ, xg, yg) 
    Sb =zeros(length(μ))
    Sbb =zeros(length(μ))
    for j in 1:length(μ)
        @show μ[j]
        bas,att = _get_pitchfork(μ[j], xg, yg)
        sb, sbb =  basin_entropy(bas)
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end

print_args = (; yticklabelsize = 8, 
            xticklabelsize = 8, 
            ylabelsize = 14, 
            xlabelsize = 14, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10",
            xticksize = 2,
            yticksize = 2)

size_cm = (8, 7); size_pt = 28.3 .* size_cm
f = Figure(resolution = size_pt)

res = 151
μ = range(0.5,1.5, length = 28)
xg = range(-1, 1, length = res)
yg = range(-1, 1, length = res)
Sb, Sbb = compute_Sb_fig(μ, xg, yg)
ax = Axis(f[1,1], ylabel = L"S_b", xlabel = L"$\mu$"; print_args...)
scatter!(ax, μ, Sb, markersize = 4, color = :black)

save("fig3b.png",f)
save("fig3b.svg",f, pt_per_unit = 1)
