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
    return @strdict(basins, att, xg, yg)
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
    @unpack bas = data
end

function compute_Sb_fig(a, J, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas,att = _get_datas(a[j], J, xg, yg)
        sb, sbb =  Attractors.basin_entropy(bas)
        t, sbb = Attractors.basins_fractal_test(bas)
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end

print_args = (; yticklabelsize = 40, 
            xticklabelsize = 40, 
            ylabelsize = 80, 
            xlabelsize = 80, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
cmap = ColorScheme([RGB(1,1,1), RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )

res = 500
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-3, 12, length = res)
a = range(1.3, 1.41, length = 40)
J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)

fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = L"S_{bb}", xlabel = L"A"; print_args...)
scatter!(ax, a, Sbb, markersize = 10, color = :black)
save("fig5_Asf.svg",fig)
save("fig5_Asf.png",fig)

# Inset 1. A1
bas, att = _get_datas(1.3, J, xg, yg)
fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], xlabel = L"x_n", ylabel = L"y_n"; print_args...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rastersize = 1)
save("fig5_inset1.png",fig)

# Inset 1. A1
bas, att = _get_datas(1.32, J, xg, yg)
fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], xlabel = L"x_n", ylabel = L"y_n"; print_args...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rastersize = 1)
save("fig5_inset2.png",fig)

# Inset 2. A1
bas, att = _get_datas(1.39, J, xg, yg)
fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], xlabel = L"x_n", ylabel = L"y_n"; print_args...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rastersize = 1)
save("fig5_inset3.png",fig)



