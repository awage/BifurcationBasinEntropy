using Distributed
@everywhere using DrWatson
@everywhere @quickactivate
@everywhere using Attractors
@everywhere using JLD2
using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes

#using Suppressor
@everywhere function henon_map!(dz, z, p, n)
    xn, yn = z
    A, J = p
    dz[1] = A - xn^2 - J*yn
    dz[2] = xn  
    return
end


@everywhere function _get_basins_henon(di)
    @unpack A, J, xg, yg = di
    u0 = [0., 0.6]
    ds = DeterministicIteratedMap(henon_map!, [1.0, 0.0], [A, J])
    xgg = range(-5, 5, length = 10000)
    ygg = range(-20, 20, length = 10000)
    grid = (xgg, ygg)
    mapper = AttractorsViaRecurrences(ds, grid,
            mx_chk_fnd_att = 5000,
            mx_chk_loc_att = 5000,
            sparse = true)
    grid = (xg, yg)
    basins, att = basins_of_attraction(mapper, grid; show_progress = true)
    sb, sbb =  basin_entropy(basins)
    t, sbb = basins_fractal_test(basins)
    return @strdict(basins, att, sb, sbb, t, xg, yg)
end

@everywhere function _get_datas(A, J, xg, yg)
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(A, J, xg, yg), # container
        _get_basins_henon, # function
        prefix = "basins_henon", # prefix for savename
        force = false,
        wsave_kwargs = (;compress = true)
    )
    @unpack att,sb,sbb = data
    return att,sb,sbb
end


function print_fig(a, P, Sb, Sbb, fname)
    print_args = (; yticklabelsize = 20, 
                xticklabelsize = 20, 
                ylabelsize = 40, 
                xlabelsize = 40, 
                xticklabelfont = "cmr10", 
                yticklabelfont = "cmr10")
    fig = Figure(resolution = (1000, 1200))
    ax = Axis(fig[1,1], ylabel = L"y"; print_args...) 
    ax2 = Axis(fig[2,1], ylabel = L"S_b"; print_args...)
    ax3 = Axis(fig[3,1], ylabel = L"S_{bb}", xlabel = L"A"; print_args...)
    pl = scatter!(ax, P[:,1],P[:,2], markersize = 0.7, color = :black, rasterize = 4)
    scatter!(ax2, a, Sb, markersize = 2.5, rasterize = 4)
    scatter!(ax3, a, Sbb, markersize = 2.5, rasterize = 4)
    Label(fig[1, 1, TopLeft()], "(a)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
    Label(fig[2, 1, TopLeft()], "(b)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
    Label(fig[3, 1, TopLeft()], "(c)", fontsize = 26, font = "cmr10", padding = (0, 5, 5, 20), halign = :right)
    save(string(fname,".svg"),fig)
    save(string(fname,".png"),fig)
end

function _get_sb(a,J)
    Sb = zeros(length(a))
    Sbb = zeros(length(a))
    for (k,ar) in enumerate(a)
        att,sb,sbb = _get_datas(a[k], J, nothing, nothing)
        Sb[k] = sb
        Sbb[k] = sbb
    end
    return Sb, Sbb
end

function compute_bif_diag(a,J)
    Sb = zeros(length(a))
    Sbb = zeros(length(a))
    pnt_lst = Vector{Vector{Float64}}(undef,1)
    df = DiscreteDynamicalSystem(henon_map!, [1.0, 0.0], [a[1], J], (J,z,p,n) -> nothing)
    for (k,ar) in enumerate(a)
        att,sb,sbb = _get_datas(a[k], J, nothing, nothing)
        Sb[k] = sb
        Sbb[k] = sbb
        if att ≠ 0 
            for p in att
                set_parameter!(df, [a[k], J])
                tra = trajectory(df, 20, p[2][1]; Ttr = 200000)
                for y in tra
                    v = [a[k], y[1]]
                    push!(pnt_lst, v)
                end
            end
        end
    end
    P = Dataset(pnt_lst[2:end])
    return P, Sb, Sbb
end

res = 5000
# SADDLE NODE 0 to 1 state. 
xg = range(-3, 3, length = res)
yg = range(-3, 12, length = res)
a = range(1.0, 2, length = 2000)
J = 0.3

# Generate bifurcation diagram 
pmap(k -> _get_datas(a[k], J, xg, yg), 1:length(a))
P, Sb, Sbb = compute_bif_diag(a,J)
Sb , Sbb = _get_sb(a,J)
save("bif_henon_J0.3.jld2", "P", P, "a", a, "J", J, "Sb", Sb, "Sbb", Sbb)

@load "bif_henon_J0.3.jld2"
print_fig(a, P, Sb, Sbb, "fig8_bif_J0.3")

