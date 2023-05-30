# Script to produce the figure for the pitchfork 
# bifurcation of the lorenz system. Paper 
# Using basin entropy to explore bifurcations
# MIT License
# Copyright (c) 2023 Alexandre Wagemkakers
# You should have received a copy of the license with this software. 
using DrWatson
@quickactivate
using Attractors
using CairoMakie
using OrdinaryDiffEq:Vern9
using LaTeXStrings

function lorenz(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3)
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e8)
    return CoupledODEs(lorenz_rule, u0, [σ, ρ, β]; diffeq = diffeq)
end
const lorenz63 = lorenz
@inbounds function lorenz_rule(u, p, t)
    du1 = p[1]*(u[2]-u[1])
    du2 = u[1]*(p[2]-u[3]) - u[2]
    du3 = u[1]*u[2] - p[3]*u[3]
    return SVector{3}(du1, du2, du3)
end
@inbounds function lorenz_jacob(u, p, t)
        return SMatrix{3,3}(-p[1], p[2] - u[3], u[2], p[1], -1.0, u[1], 0.0, -u[1], -p[3])
end

function _get_lorenz(di)
    @unpack ρ, xg, yg  = di
    ds = lorenz(rand(3); ρ = ρ)
    xgg = ygg = zgg  = range(-10, 10, length = 50001)
    grid = (xgg, ygg, zgg)
    mapper = AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 5000,
            mx_chk_loc_att = 5000)
    grid = (xg, yg)
    basins = [ mapper([x,y,0.]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    sb, sbb =  basin_entropy(basins,  20)
    return @strdict(basins, att, sb, sbb, xg, yg)
end

function _get_datas(ρ, xg, yg; force = false) 
    filename = config -> savename("basins_lorenz", config; digits = 5)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(ρ, xg, yg), # container
        _get_lorenz; # function
        filename,
        force = force,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins,att,sb,sbb = data
    return basins, att, sb, sbb
end


print_args = (; yticklabelsize = 8, 
            xticklabelsize = 8, 
            ylabelsize = 14, 
            xlabelsize = 14, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10",
            xticksize = 2,
            yticksize = 2)

size_cm = (8, 7); size_pt = 28.3 .* size_cm
f = Figure(resolution = size_pt)

res = 401
ρ = range(0.7,1.3, length = 48)
xg = range(-2, 2, length = res)
yg = range(-2, 2, length = res)
Sb =zeros(length(ρ))
Sbb =zeros(length(ρ))
for j in 1:length(ρ)
    @show ρ[j]
    bas, att, sb, sbb = _get_datas(ρ[j], xg, yg; force = true)
    Sb[j] = sb
end

ax = Axis(f[1,1], ylabel = L"S_b", xlabel = L"$\rho$"; print_args...)
scatter!(ax, ρ, Sb, markersize = 4, color = :black)
save("fig_lorenz.png",f)

