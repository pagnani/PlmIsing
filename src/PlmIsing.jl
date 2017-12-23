module PlmIsing
using NLopt, ExtractMacro
export pairplmdca, PlmOut

include("types.jl")
include("utils.jl")
include("iter.jl")
include("pairplmdca.jl")
end
