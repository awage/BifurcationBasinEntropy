@everywhere using DrWatson
@everywhere @quickactivate
@everywhere using DynamicalSystems
@everywhere using StaticArrays


@everywhere function new_henon(x, p, n)
    return SVector{2}(p[1] - x[1]^2 - (1 - p[2])*x[2],  x[1])
end



@everywhere function _get_basins_henon(d)
    @unpack res, a, ν = d
    u0 = [0., 0.6]
    df = DiscreteDynamicalSystem(new_henon, u0, [a,ν]) 
    # df = Systems.henon(u0; a, b)
    xg = yg = range(-25, 25, length = res)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(df, grid,
            mx_chk_lost = 10, 
            mx_chk_fnd_att = 3000, 
            mx_chk_loc_att = 3000, 
            mx_chk_att = 2)

    basins, att = basins_of_attraction(mapper, grid;
        show_progress = false,
)


	mn_wd, mx_wd, Sb, Sbb, α, ch_ps = 0., 0., 0., 0., 0., 0.
    Na = length(unique(basins))
	Sb , Sbb = basin_entropy(basins)

	@show a, ν, Na, Sb, Sbb

    return @strdict(basins, att, xg, yg, Na, Sb, Sbb)
end


@everywhere function compute_basins(res, a, ν)

	data, file = produce_or_load(
	    datadir("basins"), # path
	    @dict(res, a, ν), # container
	    _get_basins_henon, # function
	    prefix = "basins_henon", # prefix for savename
	    force = true,
        digits = 5,
        wsave_kwargs = (;compress = true)
	)

    return 1;
end



res = 3000
a = range(0.0, 3.98, length = 10000)
ν = 0.01

#@suppress_err begin
#cmpt_grid_wada(0.1, 0.3)
pmap(k -> compute_basins(res, a[k], ν), 1:length(a))
#end
