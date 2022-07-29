using JLD2
using Plots
using LaTeXStrings
using DynamicalSystems

function new_henon(x, p, n)
    return SVector{2}(p[1] - x[1]^2 - (1 - p[2])*x[2],  x[1])
end

@load "results_henon_res_500.jld2"
a = range(0., 3.98, length = 1000)
ν = 0.01
pnt_lst = Vector{Vector{Float64}}(undef,1)

u0 = [0., 0.6]
df = DiscreteDynamicalSystem(new_henon, u0, [a[1],ν]) 

for (k,el) in enumerate(att)
    @show el
    if el ≠ 0 
        for p in el
            set_parameter!(df, [a[k], ν])
            tra = trajectory(df, 100000, p[2][1])
        for y in tra[99990:100000]
         v = [a[k], y[1]]
         push!(pnt_lst, v)
        end
    end
end
end


P = Dataset(pnt_lst[2:end])
plot(P[:,1],P[:,2], seriestype = :scatter, markersize = 0.2, legend = false, size=(1024,768))
savefig("diag_bif3.png")

save("diag_bif.jld2", "a", a, "P", P)

