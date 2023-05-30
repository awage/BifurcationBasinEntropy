# Script to produce the snippet of the paper 
# Using basin entropy to explore bifurcations
# MIT License
# Copyright (c) 2023 Alexandre Wagemkakers
# You should have received a copy of the license with this software. 
using Attractors, DrWatson

function henon_map!(dz, z, p, n)
    xn, yn = z;   A, J = p
    dz[1] = (A - xn^2 - J*yn); dz[2] = xn 
    return 
end

function _get_basins_henon(A, J, xg, yg)
    ds = DeterministicIteratedMap(henon_map!, [1.0, 0.0], [A, J])
    xgg = range(-5, 5, length = 10000)
    ygg = range(-20, 20, length = 10000)
    mapper = AttractorsViaRecurrences(ds, (xgg,ygg),
            mx_chk_fnd_att = 500,
            mx_chk_loc_att = 500, sparse = true)
    basins, att = basins_of_attraction(mapper, (xg, yg))
    sb, sbb =  basin_entropy(basins)
    return att, sb, sbb 
end

xg = range(-3, 3, length = 1000); yg = range(-3, 12, length = 1000)
arange = range(1., 2, length = 100); J = 0.3
Sbb = zeros(100); Sb = zeros(100)
for (k,A) in enumerate(arange)
    att, sb, sbb = _get_basins_henon(A, J, xg, yg)
    Sbb[k] = sbb; Sb[k] = sb
end

