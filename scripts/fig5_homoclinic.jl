using DrWatson
@quickactivate
using DynamicalSystems
using Attractors
using StaticArrays
using JLD2
using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes
using OrdinaryDiffEq:Vern9
function bt_form!(dz, z, p, t)
    x, y = z
    λ1 = p
    dz[1] = 2*y
    dz[2] = 2*x - 3*x^2 -y*(x^3 - x^2 + y^2 -λ1) 
    return
end


function _get_bt_form(di)
    @unpack λ1, xg, yg = di
    ds = ContinuousDynamicalSystem(bt_form!, [1.0, 0.0], λ1, (J,z,p,n) -> nothing)
    xgg = range(-20, 20, length = 5000)
    ygg = range(-20, 20, length = 5000)
    grid = (xgg, ygg)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 3000, sparse = true, 
            mx_chk_loc_att = 3000,
            mx_chk_att = 2, diffeq)
    grid = (xg, yg)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  Attractors.basin_entropy(basins)
    t, sbb = Attractors.basins_fractal_test(basins)
    return @strdict(basins, att, sb, sbb, t, xg, yg)
end

function _get_datas(λ1, xg, yg)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(λ1, xg, yg), # container
        _get_bt_form, # function
        prefix = "basins_homocl", # prefix for savename
        force = false,
        digits = 5,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins, sb, sbb = data
    return basins, sb, sbb
end

function compute_Sb_fig(a, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas, sb, sbb = _get_datas(a[j], xg, yg)
        bas2 = -1 .*ones(Int16, size(bas))
        bas2[bas .== 1] .= 1
        bas2[bas .== 2] .= 1
        bas2[bas .== 3] .= 1
        sb, sbb =  Attractors.basin_entropy(bas2)
        t, sbb = Attractors.basins_fractal_test(bas2)
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

f = Figure(resolution = (1600, 600))
gb = f[1,2] = GridLayout()
gb1 = gb[1,1] = GridLayout()
gb2 = gb[1,2] = GridLayout()
ga = f[1,1] = GridLayout()


# Dummy figure for panel (a) 
ax0 = Axis(ga[1,1]; print_args...)
scatter!(ax0, [0,0], [1,1])

res = 1501
# SADDLE NODE 0 to 1 state. 
xg = range(-2, 2, length = res)
yg = range(-2, 2, length = res)
a = range(-0.2, 0.1, length = 20)
Sb, Sbb = compute_Sb_fig(a, xg, yg)

ax = Axis(gb1[1,1], ylabel = L"S_{b}", xlabel = L"c"; print_args...)
scatter!(ax, a, Sb, markersize = 10, color = :black)

bas,s,ss = _get_datas(-0.2, xg, yg)
        bas2 = -1 .*ones(Int16, size(bas))
        bas2[bas .== 1] .= 1
        bas2[bas .== 2] .= 1
        bas2[bas .== 3] .= 1
ax = Axis(gb2[1,1]; ylabel = L"y", xlabel = L"x", print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas2; colormap = cmap, rasterize = 4)


bas,s,ss = _get_datas(0.1, xg, yg)
        bas2 = -1 .*ones(Int16, size(bas))
        bas2[bas .== 1] .= 1
        bas2[bas .== 2] .= 1
        bas2[bas .== 3] .= 1
ax = Axis(gb2[2,1]; ylabel = L"y", xlabel = L"x", print_args...)
heatmap!(ax, xg, yg, bas2; colormap = cmap, rasterize = 4)

Label(ga[1, 1, TopLeft()], "(a)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gb[1, 1, TopLeft()], "(b)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
colsize!(f.layout, 1, Auto(0.7))
colsize!(gb, 2, Auto(0.5))
save("fig5.png",f)
save("fig5.svg",f)