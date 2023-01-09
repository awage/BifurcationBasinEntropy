using DynamicalSystems
using Attractors
using OrdinaryDiffEq:Vern9
using StaticArrays
using JLD2

function hopf!(dz, z, p, t)
    x, y = z
    dz[1] = -(p-x)*x 
    dz[2] = -y 
    return
end

# dummy function to keep the initializator happy
function henon_map_J(J, z0, p, n)
     return
end

function _get_hopf(res, μ)
    ds = ContinuousDynamicalSystem(hopf!, [1.0, 0.0], μ, henon_map_J)
    xg = range(-1, 1, length = res)
    yg = range(-1, 1, length = res)
    grid = (xg, yg)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = AttractorsViaRecurrences(ds, grid; Δt = 0.1,
            mx_chk_fnd_att = 100,
            mx_chk_loc_att = 100,
            mx_chk_att = 2, diffeq)
    # basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end


res = 501
a = range(-0.4,0.4, length = 20)
Sb =zeros(length(a))
Sbb =zeros(length(a))
att_m = []
for j in 1:length(a)
    @show a[j]
    bas,att = _get_hopf(res,a[j])
    sb, sbb =  basin_entropy(bas)
    Sb[j] = sb
    Sbb[j] = sbb
    push!(att_m, att)
end
plot(a, Sb)
# save("bif_henon_500.jld2", "Sb", Sb, "Sbb", Sbb, "att_m", att_m) 

