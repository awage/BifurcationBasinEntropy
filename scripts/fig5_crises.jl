# Script to produce the Figure 5 of the paper 
# Using basin entropy to explore bifurcations
# MIT License
# Copyright (c) 2023 Alexandre Wagemkakers
# You should have received a copy of the license with this software. 
using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes

#using Suppressor
function henon_map!(dz, z, p, n)
    xn, yn = z
    A, J = p
    dz[1] = A - xn^2 - J*yn
    dz[2] = xn  
    return
end


function _get_basins_henon(di)
    @unpack A, J, xg, yg = di 
    u0 = [0., 0.6]
    ds = DeterministicIteratedMap(henon_map!, [1.0, 0.0], [A, J])
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000,
            mx_chk_att = 2)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  basin_entropy(basins)
    return @strdict(basins, att, sb, sbb, xg, yg)
end

function _get_datas(A, J, xg, yg)
    filename = config -> savename("basins_henon", config; digits = 5)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(A, J, xg, yg), # container
        _get_basins_henon; # function
        filename,
        force = false,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins,att,sb,sbb = data
    return basins,att,sb,sbb
end

function compute_Sb_fig(a, J, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas,att,sb,sbb = _get_datas(a[j], J, xg, yg)
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


res = 1500
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res); yg = range(-20, 3, length = res)
a = range(1.85, 1.92, length = 40); J = 0.05
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)

ax = Axis(gb1[1,1], ylabel = L"S_b", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb, markersize = 4, color = :black)

# Inset 1. A1
bas, att, sb, sbb = _get_datas(1.8874, J, xg, yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
ax = Axis(gb2[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args1...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)

# Inset 2. A1
bas, att, sb, sbb = _get_datas(1.95, J, xg, yg)
ax = Axis(gb2[2,1]; ylabel = L"y_n", xlabel = L"x_n", print_args1...)
cmap = ColorScheme([RGB(230/255,230/255,230/255)])
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)


colsize!(gb, 2, Auto(0.35))
rowgap!(gb2, 2)
colgap!(gb, 2)
save("fig5b.png",f)
save("fig5b.svg",f, pt_per_unit = 1)

