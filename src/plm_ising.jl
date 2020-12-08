using ArgParse
using DelimitedFiles
using Distributed
@everywhere using Pkg   
@everywhere Pkg.activate(joinpath(dirname(@__FILE__),".."))
@everywhere using PlmIsing
function main()
    validmethods = ["LD_MMA", "LD_SLSQP", "LD_LBFGS", "LD_TNEWTON_PRECOND", "LD_TNEWTON_PRECOND_RESTART", "LD_TNEWTON", "LD_VAR2", "LD_VAR1"]
    s = ArgParseSettings("Ising PlmDCA inference algorithm.")
    @add_arg_table s begin
        "--lambdaJ"
        arg_type = Float64
        default = 0.01
        help = """coupling J regularization"""
        "--lambdaH"
        arg_type = Float64
        default=0.01
        help = """coupling H regularization"""
        "--epsconv", "-e"
        arg_type = Float64
        default = 1e-20
        help = """convergence criterion for the algorithm"""
        "--maxit"
        arg_type = Int
        default = 1000
        help = """maximum number of iterations"""
        "--method"
        arg_type = String
        default = "LD_LBFGS"
        range_tester = (s -> s ∈ validmethods)
        help ="Valid options are: $(join(validmethods, ", ", " and "))"
        "filein"
        help = "input filename containing the N x M configurations"
        required = true
        "fileout"
        help = "output filename in compact format (N x N matrix. The diagonal is filled with fields)\n"
        required = true
    end




    parsed_args = parse_args(ARGS, s)
    lambdaJ = parsed_args["lambdaJ"]
    lambdaH = parsed_args["lambdaH"]
    epsconv = parsed_args["epsconv"]
    maxit = parsed_args["maxit"]
    method = parsed_args["method"]
    filein = parsed_args["filein"]
    fileout = parsed_args["fileout"]


    filein == fileout && error("cowardly refusing to overwrite input file = $filein")

    res = PlmIsing.isingplm(filein,
                      lambdaJ=lambdaJ,
                      lambdaH=lambdaH,
                      epsconv=epsconv,
                      maxit=maxit,
                      method=Symbol(method))

    J = copy(res.J)
    H = res.H

    for i in eachindex(H)
        J[i,i] = H[i]
    end

    writedlm(fileout,J)
end
main()
