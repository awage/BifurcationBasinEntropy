using DrWatson
@quickactivate 
using OrdinaryDiffEq:Vern9
using Random
using Attractors
using DynamicalSystems
# using CairoMakie
using LaTeXStrings

function morris_lecar!(du, u, p, t)
    V1 = -0.00; V2 = 0.15; V4 = 0.1; VCa = 1 ;  VL = -0.5; VK = -0.7; gCa = 1.2; gK = 2; gL = 0.5; τ = 3; 
    I, V3 = p 
    V, N = u
    M(x) = 0.5*(1 + tanh((x-V1)/V2))
    G(x) = 0.5*(1 + tanh((x-V3)/V4))
    du[1] = -gCa*M(V)*(V - VCa) -gK*N*(V - VK) -gL*(V-VL) + I 
    du[2] = 1/τ*(-N + G(V))
end


function _get_basins_morris_lecar(res, I)
    # @unpack V3, I, res = d
    V3 = 0.1
    p = [I, V3]
    xg = yg = range(-1,1,length = 1000)
    df = ContinuousDynamicalSystem(morris_lecar!,rand(2), p)
    diffeq = (reltol = 1e-9,  alg = Vern9())
    mapper = AttractorsViaRecurrences(df, (xg, yg);
            stop_at_Δt = true,
            mx_chk_fnd_att = 4000, 
            mx_chk_loc_att = 4000, 
            # mx_chk_att = 2,
             sparse = true, diffeq)
    xg = range(-0.5, 0.5,length=res)
    yg = range(-0.5, 0.5,length=res)
    basins, att = basins_of_attraction(mapper, (xg, yg))
    return basins, att
end



res = 101
a = range(0.12,0.24, length = 30)
Sb =zeros(length(a))
Sbb =zeros(length(a))
att_m = []
for j in 1:length(a)
    @show a[j]
    bas,att = _get_basins_morris_lecar(res,a[j])
    sb, sbb =  basin_entropy(bas)
    Sb[j] = sb
    Sbb[j] = sbb
    push!(att_m, att)
end
using Plots
plot(a, Sb)
