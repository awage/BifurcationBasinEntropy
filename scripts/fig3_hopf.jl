using DrWatson
@quickactivate
using DynamicalSystems
using Attractors
using JLD2
using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes

function hopf!(dz, z, p, t)
    x, y = z
    a,b = p
    dz[1] = (a-x^2-y^2)*x - 0.1*y 
    dz[2] = (a-x^2-y^2)*y + 0.1*x 
    return
end


function _get_hopf(di)
    @unpack μ,xg,yg  = di
    ds = DiscreteDynamicalSystem(hopf!, [1.0, 0.0], [μ, 0.9], (J,z,p,n) -> nothing)
    xgg = range(-5, 5, length = 5001)
    ygg = range(-5, 5, length = 5001)
    grid = (xgg, ygg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; 
            mx_chk_fnd_att = 10000,
            mx_chk_loc_att = 10000,
            mx_chk_att = 2)
    grid = (xg, yg)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  Attractors.basin_entropy(basins,  20)
    return @strdict(basins, att, sb, sbb, xg, yg)
end

function _get_datas(μ, xg, yg)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(μ, xg, yg), # container
        _get_hopf, # function
        prefix = "basins_hopf", # prefix for savename
        force = false,
        digits = 5,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins,att,sb,sbb = data
    return basins,att,sb,sbb
end


function compute_Sb_fig(μ, xg, yg) 
    Sb =zeros(length(μ))
    Sbb =zeros(length(μ))
    for j in 1:length(μ)
        @show μ[j]
        bas,att, sb, sbb= _get_datas(μ[j], xg, yg)
        # sb, sbb =  Attractors.basin_entropy(bas,  20)
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
gb = f[1,2] = GridLayout()
gb1 = gb[1,1] = GridLayout()
gb2 = gb[1,2] = GridLayout()
ga = f[1,1] = GridLayout()

# Dummy figure for panel (a) 
ax0 = Axis(ga[1,1]; print_args...)
scatter!(ax0, [0,0], [1,1])

res = 4500
μ = range(-1.3, -0.7, length = 40)
xg = range(-1, 1, length = res)
yg = range(-1, 1, length = res)
Sb, Sbb = compute_Sb_fig(μ, xg, yg)

ax = Axis(gb1[1,1], ylabel = L"S_b", xlabel = L"$\mu$"; print_args...)
scatter!(ax, μ, Sb, markersize = 10, color = :black)

# Inset 1.
xg = yg =  range(-1, 1, length = 200); μ = -1.5
@unpack basins = _get_hopf(@dict(μ,xg,yg))
ax1 = Axis(gb2[1,1]; ylabel = L"y_n", xlabel = L"x_n",  print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
lp1 = heatmap!(ax1, xg, yg, basins; colormap = cmap)
lp1.rasterize = 10

xg = yg =  range(-1, 1, length = 200); μ = -0.5
@unpack basins = _get_hopf(@dict(μ,xg,yg))
ax1 = Axis(gb2[2,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
lp2 = heatmap!(ax1, xg, yg, basins; colormap = cmap)
lp2.rasterize = 10

Label(ga[1, 1, TopLeft()], "(a)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gb[1, 1, TopLeft()], "(b)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
colsize!(f.layout, 1, Auto(0.7))
colsize!(gb, 2, Auto(0.5))
save("fig3.png",f)
save("fig3.pdf",f)

