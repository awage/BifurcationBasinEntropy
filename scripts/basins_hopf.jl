using DynamicalSystems
using Attractors
using StaticArrays
using JLD2

function hopf!(dz, z, p, t)
    x, y = z
    a,b = p
    dz[1] = (a-x^2-y^2)*x - 0.1*y 
    dz[2] = (a-x^2-y^2)*y + 0.1*x 
    return
end


function _get_hopf(res, μ)
    ds = DiscreteDynamicalSystem(hopf!, [1.0, 0.0], [μ, 0.9], (J,z,p,n) -> nothing)
    xg = range(-5, 5, length = 5001)
    yg = range(-5, 5, length = 5001)
    grid = (xg, yg)
    # diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 10000,
            mx_chk_loc_att = 10000,
            mx_chk_att = 2)#, diffeq)
    xg = range(-1, 1, length = res)
    yg = range(-1, 1, length = res)
    grid = (xg, yg)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end


res = 2001
a = range(-1.5, 0, length = 50)
Sb =zeros(length(a))
Sbb =zeros(length(a))
att_m = []
for j in 1:length(a)
    @show a[j]
    bas,att = _get_hopf(res,a[j])
    sb, sbb =  Attractors.basin_entropy(bas)
    Sb[j] = sb
    Sbb[j] = sbb
    push!(att_m, att)
end
plot(a, Sb)
# save("bif_henon_500.jld2", "Sb", Sb, "Sbb", Sbb, "att_m", att_m) 

