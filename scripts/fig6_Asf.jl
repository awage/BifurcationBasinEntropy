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
    t, sbb = Attractors.basins_fractal_test(basins)
    return @strdict(basins, att, sb, sbb, t, xg, yg)
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
    @unpack basins, sb, sbb = data
    return basins, sb, sbb
end

function compute_Sb_fig(a, J, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas, sb, sbb = _get_datas(a[j], J, xg, yg)
        # sb, sbb =  Attractors.basin_entropy(bas)
        # t, sbb = Attractors.basins_fractal_test(bas)
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


f = Figure(resolution = (1600, 1000))
ga = f[1,1] = GridLayout()
gb = f[1,2] = GridLayout()
gc = f[2,1] = GridLayout()
gc1 = gc[1,1] = GridLayout()
gc2 = gc[2,1] = GridLayout()
gd = f[2,2] = GridLayout()
gd1 = gd[1,1] = GridLayout()
gd2 = gd[2,1] = GridLayout()


res = 4000
xg = range(-3, 3, length = res)
yg = range(-3, 12, length = res)

a = range(1.394, 1.397, length = 40)
J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)
ax = Axis(gd1[1,1], ylabel = L"S_{b}"; print_args...)
ax2 = Axis(gd2[1,1], ylabel = L"S_{bb}", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb, markersize = 10, color = :black)
scatter!(ax2, a, Sbb, markersize = 10, color = :black)

a = range(1.305, 1.325, length = 40)
J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)
ax = Axis(gc1[1,1], ylabel = L"S_{b}"; print_args...)
ax2 = Axis(gc2[1,1], ylabel = L"S_{bb}", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb, markersize = 10, color = :black)
scatter!(ax2, a, Sbb, markersize = 10, color = :black)




res = 1000
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-3, 12, length = res)
bas,s,ss = _get_datas(1.305, J, xg, yg)
ax = Axis(ga[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)

bas,s,ss = _get_datas(1.36, J, xg, yg)
fig = Figure(resolution = (1024, 768))
ax = Axis(gb[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)

Label(ga[1, 1, TopLeft()], "(a)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gb[1, 1, TopLeft()], "(b)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gc[1, 1, TopLeft()], "(c)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gd[1, 1, TopLeft()], "(d)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
save("fig6.png",f)
save("fig6.svg",f)

