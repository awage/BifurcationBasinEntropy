using DrWatson
@quickactivate
using DynamicalSystems
using Attractors 
using JLD2
using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes

#using Suppressor
function henon_map!(dz, z, p, n)
    xn, yn = z
    A, J = p
    dz[1] = A - xn^2 - J*yn
    dz[2] = xn  
    return
end


function _get_basins_henon(di)
    @unpack A, J, xg, yg = di
    u0 = [0., 0.6]
    ds = DiscreteDynamicalSystem(henon_map!, [1.0, 0.0], [A, J], (J,z,p,n) -> nothing)
    xgg = range(-5, 5, length = 10000)
    ygg = range(-20, 20, length = 10000)
    grid = (xgg, ygg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 5000,
            mx_chk_loc_att = 5000,
            sparse = true)
    grid = (xg, yg)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  Attractors.basin_entropy(basins)
    # t, sbb = Attractors.basins_fractal_test(basins)
    return @strdict(basins, att, sb, sbb, xg, yg)
end

function _get_datas(A, J, xg, yg)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(A, J, xg, yg), # container
        _get_basins_henon, # function
        prefix = "basins_henon", # prefix for savename
        force = false,
        digits = 5,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins,att,sb,sbb = data
    return basins,att,sb,sbb
end


function compute_Sb_fig(a, J, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas,att,sb,sbb = _get_datas(a[j], J, xg, yg)
        # sb, sbb =  Attractors.basin_entropy(bas)
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end

function print_basins(xg, yg, bas, cmap, print_args, fname)
    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1,1], xlabel = L"x_n", ylabel = L"y_n"; print_args...)
    heatmap!(ax, xg, yg, bas; colormap = cmap, rastersize = 1)
    save(fname,fig)
end

print_args = (; yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 70, 
            xlabelsize = 70, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")


f = Figure(resolution = (2000, 1200))
ga = f[1,2] = GridLayout()
ga1 = ga[1,1] = GridLayout()
ga2 = ga[1,2] = GridLayout()
gb = f[2,1] = GridLayout()
gb1 = gb[1,1] = GridLayout()
gb2 = gb[1,2] = GridLayout()
gc = f[1,1] = GridLayout()

# Dummy figure for panel (a) 
ax0 = Axis(gc[1,1]; print_args...)
scatter!(ax0, [0,0], [1,1])

# Panel (c) FOR SADDLE NODE ( period 3 within another basin): 
xg = range(-3, 3, length = res)
yg = range(-100, 80, length = res)
a = range(1.6, 1.65, length = 40)
J = 0.05
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)

ax = Axis(gb1[1,1], ylabel = L"S_b", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb, markersize = 10, color = :black)
basins, att, sb, sbb = _get_datas(1.6293745,J,xg,yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0), RGB(0,0,1), RGB(1,85/255,85/255)] )
ax4 = Axis(gb2[1,1]; print_args...)
lp2 = heatmap!(ax4, xg, yg, basins; colormap = cmap)
lp2.rasterize = 4
basins, att, sb, sbb = _get_datas(1.629375,J,xg,yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255), RGB(0,0,1)] )
ax5 = Axis(gb2[2,1]; print_args...)
lp = heatmap!(ax5, xg, yg, basins; colormap = cmap)
lp.rasterize = 4
colsize!(ga, 2, Auto(0.5))
colsize!(gb, 2, Auto(0.5))
colsize!(f.layout, 1, Auto(1))


# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-20, 3, length = res)
a = range(-0.41, -0.43, length = 40)
J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)
ax = Axis(ga1[1,1], ylabel = L"S_b", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb, markersize = 10, color = :black)
basins, att, sb, sbb = _get_datas(-0.423,J,xg,yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255)] )
ax4 = Axis(ga2[1,1]; print_args...)
lp2 = heatmap!(ax4, xg, yg, basins; colormap = cmap)
lp2.rasterize = 4
basins, att, sb, sbb = _get_datas(-0.421,J,xg,yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(0,0,1), RGB(1,85/255,85/255)] )
ax5 = Axis(ga2[2,1]; print_args...)
lp = heatmap!(ax5, xg, yg, basins; colormap = cmap)
lp.rasterize = 4
colsize!(ga, 2, Auto(0.5))
colsize!(gb, 2, Auto(0.5))
colsize!(f.layout, 1, Auto(1))
save("fig1_sn_A1_tst.png",f)
save("fig1_sn_A1_tst.pdf",f)
save("fig1_sn_A1_tst.svg",f)
