using DrWatson
@quickactivate
using DynamicalSystems
using JLD2
using Attractors
using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes

#using Suppressor
function quadratic_map!(dz, z, p, n)
    xn, yn = z
    a, b = p
    dz[1] = a*xn + yn
    dz[2] = xn^2 + b 
    return
end


function _get_basins_quadratic(di)
    @unpack a, b, xg, yg = di
    u0 = [0., 0.6]
    ds = DiscreteDynamicalSystem(quadratic_map!, [1.0, 0.0], [a, b], (J,z,p,n) -> nothing)
    xgg = range(-10, 10, length = 4001)
    ygg = range(-10, 10, length = 4001)
    grid = (xgg, ygg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid;
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000, sparse = true)
    grid = (xg,yg)
    basins, att = Attractors.basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  Attractors.basin_entropy(basins)
    # t, sbb = Attractors.basins_fractal_test(basins)
    return @strdict(basins, att, sb, sbb,  xg, yg)
end

function _get_datas(a, b, xg, yg)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(a, b, xg, yg), # container
        _get_basins_quadratic, # function
        prefix = "basins_quadratic", # prefix for savename
        force = false,
        digits = 5,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins, att, sb, sbb = data
    return basins, att, sb, sbb
end

function compute_Sb_fig_b(a, b, xg, yg) 
    Sb =zeros(length(b))
    Sbb =zeros(length(b))
    for j in 1:length(b)
        @show b[j]
        bas, att, sb, sbb = _get_datas(a, b[j], xg, yg)
        # sb, sbb =  attractors.basin_entropy(bas)
        # t, sbb = attractors.basins_fractal_test(bas)
        if isnan(sbb)
                sbb = 0
        end
        Sb[j] = sb
        Sbb[j] = sbb
    end
    return Sb, Sbb
end

function compute_bif_diag_b(a,b)
    pnt_lst = Vector{Vector{Float64}}(undef,1)
    df = DiscreteDynamicalSystem(quadratic_map!, [1.0, 0.0], [a, b[1]], (J,z,p,n) -> nothing)
    for (k,ar) in enumerate(b)
        bas,att,sb,sbb = _get_datas(a, b[k], nothing, nothing)
        if att ≠ 0 
            for p in att
                set_parameter!(df, [a, b[k]])
                tra = trajectory(df, 20, p[2][1]; Ttr = 200000)
                for y in tra
                    v = [b[k], y[1]]
                    push!(pnt_lst, v)
                end
            end
        end
    end
    P = Dataset(pnt_lst[2:end])
    return P
end

function print_fig(a, P, Sb, Sbb, fname)
print_args = (; yticklabelsize = 40, 
            xticklabelsize = 40, 
            ylabelsize = 80, 
            xlabelsize = 80, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    fig = Figure(resolution = (1024, 1024))
    # ax = Axis(fig[1,1], ylabel = L"y"; print_args...) 
    ax2 = Axis(fig[1,1], ylabel = L"S_b"; print_args...)
    ax3 = Axis(fig[2,1], ylabel = L"S_{bb}", xlabel = L"b"; print_args...)
    # pl = scatter!(ax, P[:,1],P[:,2], markersize = 0.7, color = :black)
    # pl.rasterize = true
    scatter!(ax2, a, Sb, markersize = 5)
    scatter!(ax3, a, Sbb, markersize = 5)
    # lines!(ax3, [a[1], a[end]], [0.4395093, 0.4395093]; linestyle = :dash, linewidth = 4)
    save(string(fname, ".svg"),fig)
    save(string(fname,".png"),fig)
end


print_args = (; yticklabelsize = 40, 
            xticklabelsize = 40, 
            ylabelsize = 40, 
            xlabelsize = 40, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
cmap = ColorScheme([RGB(1,1,1), RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )



f = Figure(resolution = (1600, 1000))
ga = f[1,1] = GridLayout()
ga1 = ga[1,1] = GridLayout()
ga2 = ga[2,1] = GridLayout()
gb = f[1,2] = GridLayout()
gb1 = gb[1,1] = GridLayout()
gb2 = gb[2,1] = GridLayout()


res = 1000
xg = yg = range(-3, 3, length = res)
b = range(-1.32, -1.25, length = 100)
a = -0.42
Sb, Sbb = compute_Sb_fig_b(a, b, xg, yg)
Sbb[isnan.(Sbb)] .= 0 
ax2 = Axis(ga1[1,1], ylabel = L"S_b"; print_args...)
ax3 = Axis(ga2[1,1], ylabel = L"S_{bb}", xlabel = L"b"; print_args...)
scatter!(ax2, b, Sb, markersize = 5)
scatter!(ax3, b, Sbb, markersize = 5)


xg = yg = range(-3, 3, length = res)
bas,at,s,ss = _get_datas(a, -1.32, xg, yg)
ax = Axis(gb1[1,1]; xlabel = L"$x_n$", ylabel = L"$y_n$", print_args...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)
lines!(ax, [-3, 3], [-1.32, -1.32])
lines!(ax, [0,0], [-3, 3])

bas,at,s,ss = _get_datas(a, -1.275, xg, yg)
ax = Axis(gb2[1,1]; xlabel = L"$x_n$", ylabel = L"$y_n$", print_args...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)
lines!(ax, [-3, 3], [-1.275, -1.275])
lines!(ax, [0, 0], [-3, 3])

Label(ga[1, 1, TopLeft()], "(a)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
Label(gb[1, 1, TopLeft()], "(b)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
colsize!(f.layout, 2, Auto(0.5))
save("fig7.png",f)
save("fig7.svg",f)

