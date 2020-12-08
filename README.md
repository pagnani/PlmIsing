PseudoLikelihood Maximization for Ising Model
=============================================

| **Documentation**                       | **Build Status**                                                                                | **Coverage** |
|:---------------------------------------:|:-----------------------------------------------------------------------------------------------:|:------------:|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://pagnani.github.io/PlmIsing/dev) | [![Build Status](https://travis-ci.com/pagnani/PlmIsing.svg?branch=master)](https://travis-ci.com/pagnani/PlmIsing) [![Build status](https://ci.appveyor.com/api/github/webhook?id=cyq9kwo8w2bpdgys)](https://ci.appveyor.com/project/pagnani/plmising) | [![codecov](https://codecov.io/gh/pagnani/PlmIsing/branch/master/graph/badge.svg)](https://codecov.io/gh/pagnani/PlmIsing) |

The package compute the standard pseudolikelihood optimization on pairwise spin systems (Ising). 
The method exported is ``out=isingplmdca("nomefile",kwds...)`` analyzes data stored in `nomefile` file. Data should be a `N x M` formatted file of `M` configurations of `N` spins.

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
The package is not yet registered. To use it, clone the repository locally:
```
> git clone https://github.com/pagnani/PlmIsing.git
```
Then `cd` into the `PlmIsing` directory, and  from the julia REPL do a 
```
(PlmIsing) pkg> activate .
(PlmIsing) pkg> resolve
```
If you want to use it in parallel start julia with `julia -p nprocs` where `nprocs` is the number of process 
you want to use. Then from julia REPL
```
julia> @everywhere using Pkg
julia> @everywhere Pkg.activate(".")
julia> @everywhere using PlmIsing
```

Script
------

The inference can be also run with the `plm_ising.jl` script. It requires the installation of:

* [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl):
```
julia> Pkg.add("ArgParse")
```

To use it, just run from the shell:

```
$ julia PATH-TO-PACKAGE/src/plm_ising.jl infile outfile
```
where `infile` is a file containing the Ising spin configurations (see above). Note that the script can be run in parallel just with
```
$ julia -p nproc PATH-TO-PACKAGE/src/plm_ising.jl infile outfile
```
where `nproc` is the (integer) number of cores.

Versions
--------

Version 0.1.0 works with julia v0.6
Version 0.2.0 works with julia v1.0, v1.1, v1.2
Version 0.3.0 works with julia from v1.0 ... v1.5.*
