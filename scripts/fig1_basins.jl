using DrWatson
@quickactivate
#using DynamicalSystems
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
    ds = DeterministicIteratedMap(henon_map!, [1.0, 0.0], [A, J])
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
    filename = config -> savename("basins_henon", config; digits = 5)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(A, J, xg, yg), # container
        _get_basins_henon; # function
        filename,
        force = false,
        wsave_kwargs = Dict(:compress => true)
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
    heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 1)
    save(fname,fig)
end

print_args1 = (; yticklabelsize = 8,
            xticklabelsize = 8,
            ylabelsize = 14,
            xlabelsize = 14,
            xticklabelfont = "cmr10",
            yticklabelfont = "cmr10",
            xticksize = 2,
            yticksize = 2)
print_args = (; yticklabelsize = 6,
            xticklabelsize = 6,
            ylabelsize = 8,
            xlabelsize = 8,
            xticklabelfont = "cmr10",
            xticksize = 2,
            yticksize = 2,
            yticklabelfont = "cmr10")

size_cm = (12, 7); size_pt = 28.3 .* size_cm
f = Figure(resolution = size_pt)
gb = f[1,1] = GridLayout()
gb1 = gb[1,1] = GridLayout()
gb2 = gb[1,2] = GridLayout()
g = Figure(resolution = size_pt)
gc = g[1,1] = GridLayout()
gc1 = gc[1,1] = GridLayout()
gc2 = gc[1,2] = GridLayout()

# Panel (c) FOR SADDLE NODE ( period 3 within another basin):
res = 1500
xg = range(-3, 3, length = res); yg = range(-100, 80, length = res)
a = range(1.6, 1.65, length = 40); J = 0.05
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)
ax = Axis(gc1[1,1], ylabel = L"S_b", xlabel = L"A"; print_args1...)
scatter!(ax, a, Sb, markersize = 4, color = :black)

xg = range(-2.5, -1, length = 200); yg = range(-20, 20, length = 200); A = 1.6293745;
@unpack basins = _get_basins_henon(@dict(A,J,xg,yg))
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0), RGB(0,0,1), RGB(1,85/255,85/255)] )
ax4 = Axis(gc2[2,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
heatmap!(ax4, xg, yg, basins; colormap = cmap, rasterize = 4)

xg = range(-2.5, -1, length = 200); yg = range(-20, 20, length = 200); A = 1.629375;
@unpack basins = _get_basins_henon(@dict(A,J,xg,yg))
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255), RGB(0,0,1)] )
ax5 = Axis(gc2[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
heatmap!(ax5, xg, yg, basins; colormap = cmap, rasterize = 4)

colsize!(gc, 2, Auto(0.35))
rowgap!(gc2, 5)
colgap!(gc, 5)
save("fig1c.png",g)
save("fig1c.pdf",g)
save("fig1c.svg",g, pt_per_unit = 1.)

# SADDLE NODE 0 to 1 state.
xg = range(-3, 3, length = res); yg = range(-20, 3, length = res)
a = range(-0.41, -0.43, length = 40); J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)
ax = Axis(gb1[1,1], ylabel = L"S_b", xlabel = L"A"; print_args1...)
scatter!(ax, a, Sb, markersize = 4, color = :black)

basins, att, sb, sbb = _get_datas(-0.423,J,xg,yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255)] )
ax4 = Axis(gb2[2,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
lp2 = heatmap!(ax4, xg, yg, basins; colormap = cmap, rasterize = 4)

basins, att, sb, sbb = _get_datas(-0.421,J,xg,yg)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(0,0,1), RGB(1,85/255,85/255)] )
ax5 = Axis(gb2[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
lp = heatmap!(ax5, xg, yg, basins; colormap = cmap, rasterize = 4)

colsize!(gb, 2, Auto(0.35))
rowgap!(gb2, 5)
colgap!(gb, 5)
save("fig1b.png",f)
save("fig1b.svg",f, pt_per_unit = 1.)
save("fig1b.pdf",f, pt_per_unit = 1.)
