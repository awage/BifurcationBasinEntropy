using DrWatson
@quickactivate
using DynamicalSystems
using StaticArrays
using JLD2
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
    ds = DiscreteDynamicalSystem(henon_map!, [1.0, 0.0], [A, J], (J,z,p,n) -> nothing)
    grid = (xg, yg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000,
            mx_chk_att = 2)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  Attractors.basin_entropy(basins)
    return @strdict(basins, att, sb, sbb, xg, yg)
end

function _get_datas(A, J, xg, yg)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(A, J, xg, yg), # container
        _get_basins_henon, # function
        prefix = "basins_henon", # prefix for savename
        force = false,
        digits = 5,
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
        # sb, sbb =  Attractors.basin_entropy(bas)
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end

print_args = (; yticklabelsize = 20, 
            xticklabelsize = 20, 
            ylabelsize = 40, 
            xlabelsize = 40, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
cmap = ColorScheme([RGB(1,1,1), RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )


f = Figure(resolution = (1600, 600))
gb = f[1,2] = GridLayout()
gb1 = gb[1,1] = GridLayout()
gb2 = gb[1,2] = GridLayout()
ga = f[1,1] = GridLayout()

# Dummy figure for panel (a) 
ax0 = Axis(ga[1,1]; print_args...)
scatter!(ax0, [0,0], [1,1])

res = 1500
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-20, 3, length = res)
a = range(1.85, 1.92, length = 40)
J = 0.05
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)

ax = Axis(gb1[1,1], ylabel = L"S_b", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb, markersize = 10, color = :black)

# Inset 1. A1
bas, att, sb, sbb = _get_datas(1.8874, J, xg, yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
ax = Axis(gb2[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)

# Inset 2. A1
bas, att, sb, sbb = _get_datas(1.95, J, xg, yg)
ax = Axis(gb2[2,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255)])
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)


Label(ga[1, 1, TopLeft()], "(a)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gb[1, 1, TopLeft()], "(b)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
colsize!(f.layout, 1, Auto(0.7))
colsize!(gb, 2, Auto(0.5))
save("fig4.png",f)
save("fig4.svg",f)

