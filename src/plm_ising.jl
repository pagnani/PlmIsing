#!/usr/bin/env /Users/pagnani/julia-0.6/julia
include(joinpath(dirname(@__FILE__), "../src/PlmIsing.jl"))
using .PlmIsing
function main()
    if length(ARGS) == 2     
        ifile = ARGS[1]
        v = split(ifile,"=")
        v[1] == "--in" && (filein = v[2])
        isfile(filein) || error("non existent file ", filein)
        
        ofile = ARGS[2]
        v = split(ofile,"=") 
        v[1] == "--ou" && (fileou = v[2])

        filein == fileou && error("cowardly refusing to overwrite filein = $filein")
        res = PlmIsing.isingplmdca(filein)
        J = copy(res.J)
        H = res.H

        for i in eachindex(H)
            J[i,i] = H[i]
        end
        writedlm(fileou,J)
    else
        println("Usage: julia plm_ising.jl --in=infile --ou=outfile")
        println("Output format: N x N asymmetric matrix. Fields are stored in the diagonal.")
    end
end
main()
# function main()
#     s = ArgParseSettings("Match paralog species based on co-evolution signals.")
#     @add_arg_table s begin
#         "--lambdaJ"
#         arg_type = Float64
#         help = """coupling J regularization"""
#         "--lambdaH"
#         arg_type = Real
#         help = """coupling H regularization"""
#         "--epsconv", "-e"
#         arg_type = Real
#         help = """convergence criterion for the algorithm"""
#         "--maxit"
#         arg_type = Int
#         help = """maximum number of iterations"""
#         "--method"
#         arg_type = Symbol
#         default = :LD_LBFGS
#         help ="""Default method -> :LD_LBFGS. 
#                  Valid options  -> :LD_MMA, :LD_SLSQP, :LD_LBFGS, :LD_TNEWTON_PRECOND, :LD_TNEWTON_PRECOND_RESTART, :LD_TNEWTON, :LD_VAR2, :LD_VAR1"""         
#         "filein"
#         help = "input filename containing the N x M configurations"
#         required = true
#         "fileout"
#         help = "output filename in compact format (N x N matrix. The diagonal is filled with fields)\n"
#         required = true        
#     end
    
#     parsed_args = parse_args(ARGS, s)
#     lambdaJ = parsed_args["lambdaJ"]
#     lambdaH = parsed_args["lambdaH"]
#     epsconv = parsed_args["epsconv"]
#     maxit = parsed_args["maxit"]
#     method = parsed_args["method"]
#     res = isingplmdca(filename::AbstractString,
#                 lambdaJ::Real=0.01,
#                 lambdaH::Real=0.01,
#                 epsconv::Real=1.0e-5,
#                 maxit::Int=1000,
#                 method::Symbol=:LD_LBFGS)

#     J = copy(res.J)
#     H = res.H

#     for i in eachindex(H)
#         J[i,i] = H[i]
#     end

#     writedlm(fileout,J)    
# end
# main()
