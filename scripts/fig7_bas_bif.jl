using DrWatson
@quickactivate
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
    ds = DeterministicIteratedMap(quadratic_map!, [1.0, 0.0], [a, b])
    xgg = range(-10, 10, length = 4001)
    ygg = range(-10, 10, length = 4001)
    grid = (xgg, ygg)
    mapper = AttractorsViaRecurrences(ds, grid;
            mx_chk_fnd_att = 3000,
            mx_chk_loc_att = 3000, sparse = true)
    grid = (xg,yg)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  basin_entropy(basins)
    return @strdict(basins, att, sb, sbb,  xg, yg)
end

function _get_datas(a, b, xg, yg)
    filename = config -> savename("basins_quadratic", config; digits = 5)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(a, b, xg, yg), # container
        _get_basins_quadratic; # function
        filename, 
        force = false,
        wsave_kwargs = (;compress = true)
    )
    @unpack basins, att, sb, sbb = data
    return basins, att, sb, sbb
end

function compute_Sb_fig_b(a, b, xg, yg) 
    Sb =zeros(length(b))
    Sbb =zeros(length(b))
    for j in 1:length(b)
        bas, att, sb, sbb = _get_datas(a, b[j], xg, yg)
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
    ax2 = Axis(fig[1,1], ylabel = L"S_b"; print_args...)
    ax3 = Axis(fig[2,1], ylabel = L"S_{bb}", xlabel = L"b"; print_args...)
    scatter!(ax2, a, Sb, markersize = 5)
    scatter!(ax3, a, Sbb, markersize = 5)
    save(string(fname, ".svg"),fig)
    save(string(fname,".png"),fig)
end


print_args = (; yticklabelsize = 8, 
            xticklabelsize = 8, 
            ylabelsize = 14, 
            xlabelsize = 14, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10",
            xticksize = 2,
            yticksize = 2)
print_args1 = (; yticklabelsize = 6, 
            xticklabelsize = 6, 
            ylabelsize = 8, 
            xlabelsize = 8, 
            xticklabelfont = "cmr10", 
            xticksize = 2,
            yticksize = 2,
            yticklabelfont = "cmr10")

size_cm = (12, 8); size_pt = 28.3 .* size_cm
f = Figure(resolution = size_pt)
ga = f[1,1] = GridLayout()
gb = f[1,2] = GridLayout()


res = 1000
xg = yg = range(-3, 3, length = res)
b = range(-1.32, -1.25, length = 100)
a = -0.42
Sb, Sbb = compute_Sb_fig_b(a, b, xg, yg)
Sbb[isnan.(Sbb)] .= 0 
ax2 = Axis(ga[1,1], ylabel = L"S_b"; xticklabelsvisible = false, print_args...)
ax3 = Axis(ga[2,1], ylabel = L"S_{bb}", xlabel = L"b"; print_args...)
scatter!(ax2, b, Sb, markersize = 2, color = :black)
scatter!(ax3, b, Sbb, markersize = 2, color = :black)


xg = yg = range(-3, 3, length = res)
bas,at,s,ss = _get_datas(a, -1.32, xg, yg)
ax = Axis(gb[1,1]; xlabel = L"$x_n$", ylabel = L"$y_n$", print_args1...)
cmap = ColorScheme([RGB(230/255,230/255,230/255), RGB(1,0,0),  RGB(1,85/255,85/255)] )
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)
lines!(ax, [-3, 3], [-1.32, -1.32])
lines!(ax, [0,0], [-3, 3])

bas,at,s,ss = _get_datas(a, -1.275, xg, yg)
ax = Axis(gb[2,1]; xlabel = L"$x_n$", ylabel = L"$y_n$", print_args1...)
heatmap!(ax, xg, yg, bas; colormap = cmap, rasterize = 4)
lines!(ax, [-3, 3], [-1.275, -1.275])
lines!(ax, [0, 0], [-3, 3])

colsize!(f.layout, 2, Auto(0.4))
rowgap!(ga, 5)
rowgap!(gb, 5)
colgap!(f.layout, 2)
save("fig7.png",f)
save("fig7.svg",f, pt_per_unit = 1)

