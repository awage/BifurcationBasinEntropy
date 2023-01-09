using DynamicalSystems
using StaticArrays
using JLD2

#using Suppressor
function henon_map!(dz, z, p, n)
    xn, yn = z
    A, J = p
    dz[1] = A - xn^2 - J*yn
    dz[2] = xn  
    return
end

# dummy function to keep the initializator happy
function henon_map_J(J, z0, p, n)
     return
end


function _get_basins_henon(A, J, xg, yg)
    u0 = [0., 0.6]
    ds = DiscreteDynamicalSystem(henon_map!, [1.0, 0.0], [A, J], henon_map_J)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000,
            mx_chk_att = 2)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end


res = 1500
# a = 2.065:0.05:2.8
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-4, 13, length = res)
a = range(-0.5, 1.2, length = 100)
J = 0.3
# Asf and Aff transitions
# xg = range(-3, 3, length = res)
# yg = range(-4, 13, length = res)
# a = range(1.3, 1.43, length = 100)
# J = 0.3
# FOR SADDLE NODE (within another basin): 
# xg = range(-3, 3, length = res)
# yg = range(-100, 80, length = res)
# a = range(1.62935, 1.62940, length = 20)
# J = 0.05

# bas,att = _get_basins_henon(res,a,ν)
Sb =zeros(length(a))
Sbb =zeros(length(a))
att_m = []
# Threads.@threads for j in 1:length(a)
for j in 1:length(a)
    @show a[j]
    bas,att = _get_basins_henon(a[j],J, xg, yg)
    sb, sbb =  basin_entropy(bas)
    Sb[j] = sb
    Sbb[j] = sbb
    push!(att_m, att)
end
plot(a, Sb)
# save("bif_henon_500.jld2", "Sb", Sb, "Sbb", Sbb, "att_m", att_m) 

