 using DynamicalSystems
 using StaticArrays
using Attractors
using Plots

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





 function _get_basins_henon(res, A ,J)
    u0 = [0., 0.6]
    ds = DiscreteDynamicalSystem(henon_map!, [1.0, 0.0], [A, J], henon_map_J)
    yg = range(-22, 2, length = res)
    xg = range(-3, 3, length = res)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 10000,
            mx_chk_loc_att = 10000,
            mx_chk_att = 2)

    basins, att = basins_of_attraction(mapper, grid)
    return basins, att
end


res = 800
#a = range(0.3, 1.4, length = 10000)
J = 0.3
A = 1.33

bas, att= _get_basins_henon(res, A ,J)
heatmap(bas')
