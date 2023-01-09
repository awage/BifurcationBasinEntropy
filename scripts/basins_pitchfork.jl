using DynamicalSystems
using Attractors
using OrdinaryDiffEq:Vern9
using StaticArrays
using JLD2

function pitchfork!(dz, z, p, t)
    x, y = z
    dz[1] = p*x - x^3 
    dz[2] = 0.1*y 
    return
end

# dummy function to keep the initializator happy
function henon_map_J(J, z0, p, n)
     return
end

function _get_pitchfork(res, μ)
    # ds = ContinuousDynamicalSystem(pitchfork!, [1.0, 0.0], μ, henon_map_J)
    ds = DiscreteDynamicalSystem(pitchfork!, [1.0, 0.0], μ, henon_map_J)
    xg = range(-2, 2, length = res)
    yg = range(-1, 1, length = res)
    grid = (xg, yg)
    # diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 10000,
            mx_chk_loc_att = 10000,
            mx_chk_att = 2)#, diffeq)
    xg = range(-1, 1, length = res)
    yg = range(-1, 1, length = res)
    grid = (xg, yg)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end


res = 501
a = range(0,1.5, length = 40)
Sb =zeros(length(a))
Sbb =zeros(length(a))
att_m = []
for j in 1:length(a)
    @show a[j]
    bas,att = _get_pitchfork(res,a[j])
    sb, sbb =  basin_entropy(bas)
    Sb[j] = sb
    Sbb[j] = sbb
    push!(att_m, att)
end
plot(a, Sb)
# save("bif_henon_500.jld2", "Sb", Sb, "Sbb", Sbb, "att_m", att_m) 

