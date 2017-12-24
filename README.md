PseudoLikelihood Maximization for Ising Model
=============================================

``out=isingplmdca("nomefile",kwds...)`` compute standard pseudolikelihood on data stored in `nomefile` file. Data should be a `N x M` formatted file of `M` configurations of `N` spins.

Returns ``out`` of type ``PlmOut`` with 3 fields:
1. `out.pslike` the pseudolikelihood (`Vector{Float64}`)
2. `out.J` the coupling matrix
3. `out.H` the fields



Overview
-------

The code uses
[NLopt](https://github.com/JuliaOpt/NLopt.jl) which provides a Julia
interfaces to the free/open-source [NLopt
library](http://ab-initio.mit.edu/wiki/index.php/NLopt). The program
can be run on multiple cores previous ``addprocs(nprocs)`` where
``nprocs`` should be some integer number `np` lower or equal to your
(physical) number of cores.

Install
-------
It requires the installation of:
   
* [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl).
```
julia> Pkg.add("NLopt")(https://github.com/JuliaOpt/NLopt.jl)
```
* [ExtractMacro.jl]
```
julia> Pkg.add("ExtractMacro")(https://github.com/carlobaldassi/ExtractMacro.jl)
```
* [ExtractMacro.jl]
```

To be continued ...