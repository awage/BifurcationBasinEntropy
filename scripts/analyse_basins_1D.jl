using DrWatson
@quickactivate
using DynamicalSystems
using StaticArrays




function _get_basins_henon(d)
    @unpack res, a, b = d
#     u0 = [0., 0.6]
#     df = Systems.henon(u0; a, b)
     xg = yg = range(-2, 2, length = res)
#     grid = (xg, yg)
#     mapper = AttractorsViaRecurrences(df, grid,
#             mx_chk_fnd_att = 1000,
#             mx_chk_loc_att = 1000,
#             mx_chk_att = 2)
#
#     basins, att = basins_of_attraction(mapper, grid;
#         show_progress = false,
# )
#
#
# 	mn_wd, mx_wd, Sb, Sbb, α, ch_ps = 0., 0., 0., 0., 0., 0.
#     Na = length(unique(basins))
# 	Sb , Sbb = basin_entropy(basins)

    Na = 1; Sb = 0; Sbb = 0; basins = att = 0
    @show a, b, Na, Sb, Sbb

    return @strdict(basins, att, xg, yg, Na, Sb, Sbb)
end


function compute_basins(res, a, b)

	data, file = produce_or_load(
	    datadir("basins"), # path
	    @dict(res, a, b), # container
	    _get_basins_henon, # function
	    prefix = "basins_henon", # prefix for savename
	    force = false,
        digits = 5,
        wsave_kwargs = (;compress = true)
	)

    return data;
end



res = 500
a = range(0.3, 1.4, length = 10000)
b = 0.3



Na_mat = zeros(length(a))
Sb_mat = zeros(length(a))
Sbb_mat = zeros(length(a))


for k in 1:length(a)
	data = compute_basins(res, a[k], b)
	@unpack Na, Sb, Sbb = data
	Na_mat[k] = Na
	Sb_mat[k] = Sb
	Sbb_mat[k] = Sbb
end



save("results_henon_res_500.jld2", "Na", Na_mat, "Sb", Sb_mat, "Sbb", Sbb_mat, "a", a, "b", b)
