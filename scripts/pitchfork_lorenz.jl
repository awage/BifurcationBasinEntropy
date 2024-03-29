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

function lorenz(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3)
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e8)
    return CoupledODEs(lorenz_rule, u0, [σ, ρ, β]; diffeq = diffeq)
end
const lorenz63 = lorenz
@inbounds function lorenz_rule(u, p, t)
    du1 = p[1]*(u[2]-u[1])
    du2 = u[1]*(p[2]-u[3]) - u[2]
    du3 = u[1]*u[2] - p[3]*u[3]
    return SVector{3}(du1, du2, du3)
end
@inbounds function lorenz_jacob(u, p, t)
        return SMatrix{3,3}(-p[1], p[2] - u[3], u[2], p[1], -1.0, u[1], 0.0, -u[1], -p[3])
end

function _get_lorenz(di)
    @unpack ρ, xg, yg  = di
    ds = lorenz(rand(3); ρ = ρ)
    xgg = ygg = zgg  = range(-10, 10, length = 50001)
    grid = (xgg, ygg, zgg)
    mapper = AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 5000,
            mx_chk_loc_att = 5000)
    grid = (xg, yg)
    basins = [ mapper([x,y,0.]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    sb, sbb =  basin_entropy(basins,  20)
    return @strdict(basins, att, sb, sbb, xg, yg)
end

function _get_datas(ρ, xg, yg; force = false) 
    filename = config -> savename("basins_lorenz", config; digits = 5)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(ρ, xg, yg), # container
        _get_lorenz; # function
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

res = 401
ρ = range(0.7,1.3, length = 48)
xg = range(-2, 2, length = res)
yg = range(-2, 2, length = res)
Sb =zeros(length(ρ))
Sbb =zeros(length(ρ))
for j in 1:length(ρ)
    @show ρ[j]
    bas, att, sb, sbb = _get_datas(ρ[j], xg, yg; force = false)
    Sb[j] = sb
end


ax = Axis(gb1[1,1], ylabel = L"S_b", xlabel = L"$\rho$"; print_args...)
scatter!(ax, ρ, Sb, markersize = 4, color = :black)

# Inset 1. A1
bas, att, sb, sbb = _get_datas(0.7, xg, yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255)])
ax = Axis(gb2[2,1]; ylabel = L"y", xlabel = L"x", print_args1...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)

# Inset 2. A1
bas, att, sb, sbb = _get_datas(1.3, xg, yg)
ax = Axis(gb2[1,1]; ylabel = L"y", xlabel = L"x", print_args1...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)


colsize!(gb, 2, Auto(0.35))
rowgap!(gb2, 2)
colgap!(gb, 2)
save("fig_3c.png",f)
save("fig_3c.svg",f, pt_per_unit = 1)

