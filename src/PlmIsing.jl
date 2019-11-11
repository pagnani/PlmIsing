module PlmIsing
using NLopt
using ExtractMacro
using SharedArrays
using Distributed
using Printf
using DelimitedFiles

export isingplmdca, PlmOut

include("types.jl")
include("utils.jl")
include("iter.jl")
include("pairplmdca.jl")
end
