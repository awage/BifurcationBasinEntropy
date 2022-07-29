 using DynamicalSystems
 using StaticArrays


#using Suppressor


 function _get_basins_henon(res, a ,b)
    u0 = [0., 0.6]
    df = Systems.henon(u0; a, b)
    xg = yg = range(-5, 5, length = res)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(df, grid,
            mx_chk_fnd_att = 10000,
            mx_chk_loc_att = 10000,
            mx_chk_att = 2)

    basins, att = basins_of_attraction(mapper, grid)
    return basins, att
end


res = 400
#a = range(0.3, 1.4, length = 10000)
a = 0.3
b = 0.99

bas, att= _get_basins_henon(res, a ,b)
heatmap(bas')
#@suppress_err begin
#cmpt_grid_wada(0.1, 0.3)
#pmap(k -> compute_basins(res, a[k], b), 1:length(a))
#end
