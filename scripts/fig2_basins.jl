using DrWatson
@quickactivate
using DynamicalSystems
using Attractors
using JLD2
using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes

function pitchfork!(dz, z, p, t)
    x, y = z
    dz[1] = p*x - x^3 
    dz[2] = 0.1*y 
    return
end

function _get_pitchfork(μ, xg, yg)
    ds = DiscreteDynamicalSystem(pitchfork!, [1.0, 0.0], μ, (J,z,p,n) -> nothing)
    xgg = range(-2, 2, length = 5000)
    ygg = range(-1, 1, length = 5000)
    grid = (xgg, ygg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 10000,
            mx_chk_loc_att = 10000,
            mx_chk_att = 2)
    grid = (xg, yg)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    return basins,att 
end

function compute_Sb_fig(μ, xg, yg) 
    Sb =zeros(length(μ))
    Sbb =zeros(length(μ))
    for j in 1:length(μ)
        @show μ[j]
        bas,att = _get_pitchfork(μ[j], xg, yg)
        sb, sbb =  Attractors.basin_entropy(bas)
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end

print_args = (; yticklabelsize = 20, 
            xticklabelsize = 20, 
            ylabelsize = 40, 
            xlabelsize = 40, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
cmap = ColorScheme([RGB(1,1,1), RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )


f = Figure(resolution = (1600, 600))
ga = f[1,1] = GridLayout()
gb = f[1,2] = GridLayout()

# Dummy figure for panel (a) 
ax0 = Axis(ga[1,1]; print_args...)
scatter!(ax0, [0,0], [1,1])

res = 151
μ = range(0.5,1.5, length = 28)
xg = range(-1, 1, length = res)
yg = range(-1, 1, length = res)
Sb, Sbb = compute_Sb_fig(μ, xg, yg)
ax = Axis(gb[1,1], ylabel = L"S_b", xlabel = L"$\mu$"; print_args...)
scatter!(ax, μ, Sb, markersize = 10, color = :black)

Label(ga[1, 1, TopLeft()], "(a)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gb[1, 1, TopLeft()], "(b)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
 colsize!(f.layout, 1, Auto(1))
save("fig2.png",f)
save("fig2.svg",f)
