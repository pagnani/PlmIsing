var documenterSearchIndex = {"docs":
[{"location":"#","page":"PlmIsing Documentation","title":"PlmIsing Documentation","text":"    CurrentModule = PlmIsing\n    DocTestSetup = quote\n    using PlmIsing\nend","category":"page"},{"location":"#","page":"PlmIsing Documentation","title":"PlmIsing Documentation","text":"    Modules = [PlmIsing]","category":"page"},{"location":"#PlmIsing.isingplm-Tuple{Array{T,2} where T,Array{T,1} where T}","page":"PlmIsing Documentation","title":"PlmIsing.isingplm","text":"isingplm(spin::Matrix,W::Vector,kwds...)\n\nPseudo likelihood maximization of an ising pair spin interaction system, from a  N × M spin::Matrix where N is number of Ising spins (±1) and M is the number  of configurations. W is a vector of M non negative normalized (sum(W) ≈ 1) weights. Keyword args:\n\n* lamdaJ::Real=0.01; lagrange multiplier for the Js\n\n* lambdaH::Real=0.01; lagrange multiplier for the hs\n\n* epsconv::Real::1e-5; convergence paramter\n\n* maxeval::Int=1000; maximal number of iterations\n\n* verbose::Bool=true; verbosity of the output\n\n* method::Symbol=:LD_LBFGS; NLopt maximization strategy\n\nOutput: a PlmOut struct containing the Js (as a N × N matrix ), the h (as a length N vector),  and the site pseudo-likelihood.\n\nisingplm(filename::String;kwds...)\n\nRead the configurations from a filename. The file is a N × M matrix. In this case the weights are assumed to  be all 1.0/M.\n\n\n\n\n\n","category":"method"},{"location":"#PlmIsing-Documentation-1","page":"PlmIsing Documentation","title":"PlmIsing Documentation","text":"","category":"section"},{"location":"#","page":"PlmIsing Documentation","title":"PlmIsing Documentation","text":"    isingplm","category":"page"},{"location":"#PlmIsing.isingplm","page":"PlmIsing Documentation","title":"PlmIsing.isingplm","text":"isingplm(spin::Matrix,W::Vector,kwds...)\n\nPseudo likelihood maximization of an ising pair spin interaction system, from a  N × M spin::Matrix where N is number of Ising spins (±1) and M is the number  of configurations. W is a vector of M non negative normalized (sum(W) ≈ 1) weights. Keyword args:\n\n* lamdaJ::Real=0.01; lagrange multiplier for the Js\n\n* lambdaH::Real=0.01; lagrange multiplier for the hs\n\n* epsconv::Real::1e-5; convergence paramter\n\n* maxeval::Int=1000; maximal number of iterations\n\n* verbose::Bool=true; verbosity of the output\n\n* method::Symbol=:LD_LBFGS; NLopt maximization strategy\n\nOutput: a PlmOut struct containing the Js (as a N × N matrix ), the h (as a length N vector),  and the site pseudo-likelihood.\n\nisingplm(filename::String;kwds...)\n\nRead the configurations from a filename. The file is a N × M matrix. In this case the weights are assumed to  be all 1.0/M.\n\n\n\n\n\n","category":"function"}]
}