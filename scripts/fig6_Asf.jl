using DrWatson
@quickactivate
using Attractors
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
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000,
            mx_chk_att = 2)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  basin_entropy(basins)
    t, sbb = basins_fractal_test(basins)
    return @strdict(basins, att, sb, sbb, t, xg, yg)
end

function _get_datas(A, J, xg, yg)
    filename = config -> savename("basins_henon", config; digits = 5)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(A, J, xg, yg), # container
        _get_basins_henon; # function
        filename,
        force = false,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins, sb, sbb = data
    return basins, sb, sbb
end

function compute_Sb_fig(a, J, xg, yg) 
    Sb =zeros(length(a))
    Sbb =zeros(length(a))
    for j in 1:length(a)
        @show a[j]
        bas, sb, sbb = _get_datas(a[j], J, xg, yg)
        # sb, sbb =  Attractors.basin_entropy(bas)
        # t, sbb = Attractors.basins_fractal_test(bas)
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end


print_args = (; yticklabelsize = 8, 
            xticklabelsize = 8, 
            ylabelsize = 14, 
            xlabelsize = 14, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10",
            xticksize = 2,
            yticksize = 2)

size_cm = (7.5, 6); size_pt = 28.3 .* size_cm


res = 4000
xg = range(-3, 3, length = res)
yg = range(-3, 12, length = res)

a = range(1.394, 1.397, length = 40)
J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)
f = Figure(resolution = size_pt)
g = f[1,1] = GridLayout()
ax = Axis(g[1,1], ylabel = L"S_{b}"; xticklabelsvisible = false, print_args...)
ax2 = Axis(g[2,1], ylabel = L"S_{bb}", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb;  markersize = 4, color = :black)
scatter!(ax2, a, Sbb, markersize = 4, color = :black)

rowgap!(g, 5)
save("fig6d.png",f)
save("fig6d.svg",f, pt_per_unit = 1)


a = range(1.305, 1.325, length = 40)
J = 0.3
Sb, Sbb = compute_Sb_fig(a, J, xg, yg)
f = Figure(resolution = size_pt)
g = f[1,1] = GridLayout()
ax = Axis(g[1,1], ylabel = L"S_{b}"; xticklabelsvisible = false, print_args...)
ax2 = Axis(g[2,1], ylabel = L"S_{bb}", xlabel = L"A"; print_args...)
scatter!(ax, a, Sb, markersize = 4, color = :black)
scatter!(ax2, a, Sbb, markersize = 4, color = :black)

rowgap!(g, 5)
save("fig6c.png",f)
save("fig6c.svg",f, pt_per_unit = 1)



size_cm = (6, 5); size_pt = 28 .* size_cm
res = 1000
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-3, 12, length = res)
bas,s,ss = _get_datas(1.305, J, xg, yg)
f = Figure(resolution = size_pt)
ax = Axis(f[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)
save("fig6a.png",f)
save("fig6a.svg",f, pt_per_unit = 1)

bas,s,ss = _get_datas(1.36, J, xg, yg)
fig = Figure(resolution = (1024, 768))
f = Figure(resolution = size_pt)
ax = Axis(f[1,1]; ylabel = L"y_n", xlabel = L"x_n", print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)
save("fig6b.png",f)
save("fig6b.svg",f, pt_per_unit = 1)

