@everywhere using DrWatson
@everywhere @quickactivate
@everywhere using DynamicalSystems
@everywhere using StaticArrays


#using Suppressor


@everywhere function _get_basins_henon(d)
    @unpack res, a, b = d
    u0 = [0., 0.6]
    df = Systems.henon(u0; a, b)
    xg = yg = range(-5, 5, length = res)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(df, grid,
            mx_chk_fnd_att = 2000,
            mx_chk_loc_att = 2000,
            mx_chk_att = 2)

    basins, att = basins_of_attraction(mapper, grid;
        show_progress = false,
)


	mn_wd, mx_wd, Sb, Sbb, α, ch_ps = 0., 0., 0., 0., 0., 0.
    Na = length(unique(basins))
	Sb , Sbb = basin_entropy(basins)

	@show a, b, Na, Sb, Sbb

    return @strdict(basins, att, xg, yg, Na, Sb, Sbb)
end


@everywhere function compute_basins(res, a, b)

	data, file = produce_or_load(
	    datadir("basins"), # path
	    @dict(res, a, b), # container
	    _get_basins_henon, # function
	    prefix = "basins_henon", # prefix for savename
	    force = true,
        digits = 5,
        wsave_kwargs = (;compress = true)
	)

    return 1;
end



res = 500
a = range(0.3, 1.4, length = 10000)
b = 0.3

#@suppress_err begin
#cmpt_grid_wada(0.1, 0.3)
    pmap(k -> compute_basins(res, a[k], b), 1:length(a))
#end
