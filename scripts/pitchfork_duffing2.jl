# Script to produce the figure for the pitchfork 
# bifurcation of the lorenz system. Paper 
# Using basin entropy to explore bifurcations
# MIT License
# Copyright (c) 2023 Alexandre Wagemkakers
# You should have received a copy of the license with this software. 
using DrWatson
@quickactivate
using Attractors
using CairoMakie
using OrdinaryDiffEq:Vern9
using LaTeXStrings
using Colors
using ColorSchemes

function duffing(u0 = [0.1, 0.25]; ω = 2.2, f = 27.0, d = 0.2, β = 1)
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e8)
    return CoupledODEs(duffing_rule, u0, [ω, f, d, β]; diffeq = diffeq)
end
@inbounds function duffing_rule(x, p, t)
    ω, f, d, β = p
    dx1 = x[2]
    dx2 = f*cos(ω*t) + β*x[1] - x[1]^3 - d * x[2]
    return SVector(dx1, dx2)
end

function _get_duffing(di)
    @unpack β, xg, yg  = di
    ds = duffing(rand(2); β = β, f = 0.)
    xgg = ygg = range(-10, 10, length = 5001)
    grid = (xgg, ygg)
    mapper = AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 5000,
            mx_chk_loc_att = 5000)
    grid = (xg, yg)
    # basins = [ mapper([x,y]) for x in xg, y in yg]
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    # att = extract_attractors(mapper)
    sb, sbb =  basin_entropy(basins,  20)
    return @strdict(basins, att, sb, sbb, xg, yg)
end

function _get_datas(β, xg, yg; force = false) 
    filename = config -> savename("basins_duffing", config; digits = 5)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(β, xg, yg), # container
        _get_duffing; # function
        filename,
        force = force,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins,att,sb,sbb = data
    return basins, att, sb, sbb
end


print_args = (; yticklabelsize = 8, 
            xticklabelsize = 8, 
            ylabelsize = 14, 
            xlabelsize = 14, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10",
            xticksize = 2,
            yticksize = 2)
print_args1 = (; yticklabelsize = 6, 
            xticklabelsize = 6, 
            ylabelsize = 8, 
            xlabelsize = 8, 
            xticklabelfont = "cmr10", 
            xticksize = 2,
            yticksize = 2,
            yticklabelfont = "cmr10")

size_cm = (12, 7); size_pt = 28 .* size_cm
f = Figure(resolution = size_pt)
gb = f[1,1] = GridLayout()
gb1 = gb[1,1] = GridLayout()
gb2 = gb[1,2] = GridLayout()

res = 1501
β = range(-0.5, 0.5, length = 60)
xg = range(-2, 2, length = res)
yg = range(-2, 2, length = res)
Sb =zeros(length(β))
Sbb =zeros(length(β))
for j in 1:length(β)
    @show β[j]
    bas, att, sb, sbb = _get_datas(β[j], xg, yg; force = false)
    Sb[j] = sb
end


ax = Axis(gb1[1,1], ylabel = L"S_b", xlabel = L"$\beta$"; print_args...)
scatter!(ax, β, Sb, markersize = 4, color = :black)

# Inset 1. A1
bas, att, sb, sbb = _get_datas(-0.5, xg, yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255)])
ax = Axis(gb2[2,1]; ylabel = L"y", xlabel = L"x", print_args1...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)

# Inset 2. A1
bas, att, sb, sbb = _get_datas(0.008474, xg, yg)
ax = Axis(gb2[1,1]; ylabel = L"y", xlabel = L"x", print_args1...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)


colsize!(gb, 2, Auto(0.35))
rowgap!(gb2, 2)
colgap!(gb, 2)
save("fig_3c.png",f)
save("fig_3c.svg",f, pt_per_unit = 1)

