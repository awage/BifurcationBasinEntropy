using JLD2
using LaTeXStrings
using DynamicalSystems
using CairoMakie

function new_henon(x, p, n)
    return SVector{2}(p[1] - x[1]^2 - (1 - p[2])*x[2],  x[1])
end

@load "bif_henon_500.jld2"
# a = range(0., 4, length = 4000)
a = 0:0.001:4
ν = 0.01
pnt_lst = Vector{Vector{Float64}}(undef,1)

u0 = [0., 0.6]
df = DiscreteDynamicalSystem(new_henon, u0, [a[1],ν]) 

for (k,el) in enumerate(att_m)
    @show el
    if el ≠ 0 
        for p in el
            set_parameter!(df, [a[k], ν])
            tra = trajectory(df, 20, p[2][1]; Ttr = 200000)
        for y in tra
         v = [a[k], y[1]]
         push!(pnt_lst, v)
        end
    end
end
end

P = Dataset(pnt_lst[2:end])

save("diag_bif.jld2", "a", a, "P", P)

fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = "xn", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20)
ax2 = Axis(fig[2,1], ylabel = "Sb", xlabel = "a", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)

#ax = Axis(fig[1,1], ylabel = "y", yticklabelsvisible = false, xticklabelsvisible = false, bottomspinevisible = false, topspinevisible = false, leftspinevisible = false, rightspinevisible = false, xgridvisible = false, ygridvisible = false)
#ax2 = Axis(fig[2,1], ylabel = "Sb", xlabel = "F")

ind = findall(P[:,1] .≤ 4)
scatter!(ax, P[ind,1],P[ind,2], markersize = 0.7, color = :black, rasterize = 4)

ind = findall(a .≤ 4)

lines!(ax2, a[ind], vec(Sb[ind]))
#lp.rasterize = 4

save("diag_bif_henon.svg",fig)
run(`inkscape diag_bif_henon.svg --export-type=pdf`)

# plot(P[:,1],P[:,2], seriestype = :scatter, markersize = 0.2, legend = false, size=(1024,768))
# savefig("diag_bif3.png")


