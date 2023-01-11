using DynamicalSystems
using StaticArrays
using JLD2
using CairoMakie

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
    mapper = AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000,
            mx_chk_att = 2)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end

function compute_Sb_fig(a, J, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas,att = _get_basins_henon(a[j], J, xg, yg)
        sb, sbb =  basin_entropy(bas)
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end


res = 1500
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-20, 3, length = res)
a = range(-0.423, -0.421, length = 10)
J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)

fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = "Sb", xlabel = "A", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
lines!(ax, a, Sb, markersize = 0.7, color = :black)
save("fig1_sn_A1.svg",fig)
save("fig1_sn_A1.png",fig)

# FOR SADDLE NODE ( period 3 within another basin): 
# xg = range(-3, 3, length = res)
# yg = range(-100, 80, length = res)
# a = range(1.6293745, 1.629375, length = 10)
# J = 0.05



