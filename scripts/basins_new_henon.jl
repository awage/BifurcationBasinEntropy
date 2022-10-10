using DynamicalSystems
using StaticArrays
using JLD2

 function new_henon(x, p, n)
    return SVector{2}(p[1] - x[1]^2 - (1 - p[2])*x[2],  x[1])
end



function _get_basins_henon(res, a, ν)
    u0 = [0., 0.6]
    df = DiscreteDynamicalSystem(new_henon, u0, [a,ν]) 
    # df = Systems.henon(u0; a, b)
    xg = yg = range(-2.5, 2.5, length = res)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(df, grid,
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000,
            mx_chk_att = 1)
    # basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    basins, att = basins_of_attraction(mapper; show_progress = true)
    return basins,att 
end


res = 500
a = 0.:0.001:4
ν = 0.01

# bas,att = _get_basins_henon(res,a,ν)
Sb =zeros(length(a))
Sbb =zeros(length(a))
att_m = []
# Threads.@threads for j in 1:length(a)
for j in 1:length(a)
    @show a[j]
    bas,att = _get_basins_henon(res,a[j],ν)
    sb, sbb =  basin_entropy(bas)
    Sb[j] = sb
    Sbb[j] = sbb
    push!(att_m, att)
end

save("bif_henon_500.jld2", "Sb", Sb, "Sbb", Sbb, "att_m", att_m) 

