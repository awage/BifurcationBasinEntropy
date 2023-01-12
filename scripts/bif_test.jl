using Revise, Plots
using BifurcationKit, Setfield, ForwardDiff
const BK = BifurcationKit

FbpSecBif(u, p) = @. -u * (p + u * (2-5u)) * (p -.15 - u * (2+20u))
dFbpSecBif(x,p) =  ForwardDiff.jacobian( z-> FbpSecBif(z,p), x)
# bifurcation problem
prob = BifurcationProblem(FbpSecBif, [0.0], -0.2,
	# specify the continuation parameter
	(@lens _), J = dFbpSecBif, recordFromSolution = (x, p) -> x[1])

# options for Krylov-Newton
opt_newton = NewtonPar(tol = 1e-9, maxIter = 20)

# options for continuation
opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds = 0.01,
	maxSteps = 100, nev = 2, newtonOptions = opt_newton,
	# parameter interval
	pMax = 0.4, pMin = -0.5,
	# detect bifurcations with bisection method
	detectBifurcation = 3, nInversion = 4, tolBisectionEigenvalue = 1e-8, dsminBisection = 1e-9)

diagram = bifurcationdiagram(prob, PALC(),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	2,
	(args...) -> setproperties(opts_br; pMin = -1.0, pMax = .3, ds = 0.001, dsmax = 0.005, nInversion = 8, detectBifurcation = 3, dsminBisection =1e-18, maxBisectionSteps=20))

# You can plot the diagram like
plot(diagram; putspecialptlegend=false, markersize=2, plotfold=false, title = "#branches = $(size(diagram))")
