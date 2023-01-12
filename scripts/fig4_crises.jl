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


function _get_basins_henon(A, J, xg, yg)
    u0 = [0., 0.6]
    ds = DiscreteDynamicalSystem(henon_map!, [1.0, 0.0], [A, J], (J,z,p,n) -> nothing)
    grid = (xg, yg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000,
            mx_chk_att = 2)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end

function compute_Sb_fig(a, J, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas,att = _get_basins_henon(a[j], J, xg, yg)
        sb, sbb =  Attractors.basin_entropy(bas)
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

res = 1500
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-20, 3, length = res)
a = range(1.85, 1.92, length = 40)
J = 0.05
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)

fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = L"S_b", xlabel = L"A"; print_args...)
# lines!(ax, a, Sb, markersize = 0.7, color = :black)
scatter!(ax, a, Sb, markersize = 10, color = :black)
save("fig4_crisis.svg",fig)
save("fig4_crisis.png",fig)

# Inset 1. A1
bas, att = _get_basins_henon(1.8874, J, xg, yg)
fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], xlabel = L"x_n", ylabel = L"y_n"; print_args...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rastersize = 1)
save("fig4_cris_inset1.png",fig)

# Inset 2. A1
bas, att = _get_basins_henon(1.95, J, xg, yg)
fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], xlabel = L"x_n", ylabel = L"y_n"; print_args...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rastersize = 1)
save("fig4_cris_inset2.png",fig)



