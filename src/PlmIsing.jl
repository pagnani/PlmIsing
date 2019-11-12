module PlmIsing
using NLopt
using ExtractMacro
using SharedArrays
using Distributed
using Printf
using DelimitedFiles

export isingplm, PlmOut

include("types.jl")
include("utils.jl")
include("isingplm.jl")
include("pairplmdca.jl")
end #end module
