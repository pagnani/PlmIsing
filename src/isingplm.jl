"""
    isingplm(spin::Matrix,W::Vector,kwds...)
Pseudo likelihood maximization of an ising pair spin interaction system, from a 
N × M `spin::Matrix` where N is number of Ising spins (±1) and M is the number 
of configurations. `W` is a vector of M non negative normalized (`sum(W) ≈ 1`) weights.
Keyword args:

    * lamdaJ::Real=0.01; lagrange multiplier for the Js

    * lambdaH::Real=0.01; lagrange multiplier for the hs

    * epsconv::Real::1e-5; convergence paramter

    * maxeval::Int=1000; maximal number of iterations

    * verbose::Bool=true; verbosity of the output

    * method::Symbol=:LD_LBFGS; NLopt maximization strategy

Output: a `PlmOut` struct containing the Js (as a N × N matrix ), the h (as a length N vector), 
and the site pseudo-likelihood.


    isingplm(filename::String;kwds...) 

Read the configurations from a `filename`. The file is a N × M matrix. In this case the weights are assumed to 
be all 1.0/M.
"""
function isingplm(spin::Matrix,W::Vector;
                lambdaJ::Real=0.01,
                lambdaH::Real=0.01,
                epsconv::Real=1.0e-5,
                maxit::Int=1000,
                maxeval::Int=5000,
                verbose::Bool=true,
                method::Symbol=:LD_LBFGS)

    N,M = size(spin)
    length(W) == M || throw(DimensionMismatch("incompatible weigth vector length"))
    all(x->x>0,W) || throw(DomainError("vector W should normalized and with all positive elements"))
    isapprox(sum(W),1) || throw(DomainError("sum(W) ≠ 1. Consider normalizing the vector W"))
    plmalg = PlmAlg(method, verbose, epsconv, maxit, maxeval)
    plmvar = PlmVar(M, N, lambdaJ, lambdaH, spin,W)
    outJ, outh, vecps=maximizeplmdca(plmalg,plmvar)
    return PlmOut(sdata(vecps),outJ,outh)
end

isingplm(spin::Matrix; kwds...) = isingplm(spin,ones(size(spin,2))/size(spin,2);kwds...)

function isingplm(filename::AbstractString; kwds...)
    spin = readdlm(filename,Int)
    isingplm(spin; kwds...)
end

function optimfunwrapper(x::Vector, g::Vector, site, var)
    g === nothing && (g = zeros(Float64, length(x)))
    return plmsitegrad!(x, g, site,  var)
end

function maximizeplmdca(alg::PlmAlg, var::PlmVar)

    @extract var M N lambdaJ lambdaH spin
    @extract alg method verbose epsconv maxit maxeval
    vecps = SharedArray{Float64}(N)
    x0 = fill(0.0,N)
    Jscra = @distributed hcat for i in 1:N
        opt = Opt(method,N)
        ftol_abs!(opt, alg.epsconv)
        xtol_rel!(opt, alg.epsconv)
        xtol_abs!(opt, alg.epsconv)
        ftol_rel!(opt, alg.epsconv)
        maxeval!(opt, maxit)
        max_objective!(opt,(x,g)->optimfunwrapper(x,g,i,var))
        elapsedtime = @elapsed (maxf, maxJH, ret) = optimize(opt,x0)
        alg.verbose && @printf("site %d\t pll = %.4f\t time(s) = %.4f\t", i, maxf, elapsedtime)
        alg.verbose && println("exit status = $ret")
        #alg.verbose && println(maxJH[end-2], " ",maxJH[end-1]," ", maxJH[end], " ", mean(maxJH[1:2N-3]))
        vecps[i] = maxf
        maxJH
    end
    outJ,outh = unpack(Jscra)
    return outJ, outh, vecps
end

function computeH(vecJ::Vector{Float64}, site::Int, spin, a::Int,  N::Int)
    @inbounds begin 
        ctr = 0
        # begin hi
        Hi = vecJ[end] # the i-th field
        @simd for i=1:site-1
            ctr += 1
            Hi +=  vecJ[ctr] * spin[i,a]
        end
        @simd for i=site+1:N
            ctr += 1
            Hi += vecJ[ctr] * spin[i,a]
        end
    end
    return Hi
end

function L2norm(v::Vector, var::PlmVar)

    #pl = pl0 - lambdaJ |J|₂ - lambdaH |H|₂
    @extract var N lambdaJ lambdaH
    mysum = 0.0
    for i=1:N-1
        mysum += v[i]^2
    end
    return -lambdaJ*mysum - 2lambdaH*v[end]^2
end

function grad!(grad::Vector{Float64}, spin, Wa, N::Int, Hi::Float64, site::Int, a::Int)
    exponent = -2.0*Hi*spin[site,a]

    factor = 2.0 * spin[site,a]
    if exponent < 500.0
        factor = 2.0 * spin[site,a] * (exp(exponent) / (1.0 + exp(exponent)))
    end

    @inbounds begin
        ctr = 0
        @simd for i=1:site-1
            ctr += 1
            grad[ctr] += Wa*spin[i,a]*factor
        end
        @simd for i=site+1:N
            ctr += 1
            grad[ctr] += Wa*spin[i,a]*factor
        end
        grad[end] += Wa*factor
    end
end

function plmsitegrad!(vecJ::Vector{Float64}, grad::Vector{Float64}, site::Int, var::PlmVar)
    @extract var M N lambdaJ lambdaH

    spin = sdata(var.spin)
    W = sdata(var.W)
    pl = 0.0

    @inbounds @simd for i in 1:N-1
        grad[i] = -2 * lambdaJ * vecJ[i] 
    end
    grad[end] = -4 * lambdaH * vecJ[end] 
    
    @inbounds for a in 1:M
        #println("W[$a]= ",W[a], " M = $M size(spin) = $(size(spin)) site = $site $(typeof(spin))")
        Hi = computeH(vecJ, site, spin, a, N)
        exponent = -2.0*spin[site,a]*Hi
        plterm = -exponent
        if exponent < 700.0
            plterm = -log(1.0 +exp(exponent))
        end  
        pl += W[a] * plterm
        grad!(grad, spin, W[a] , N, Hi, site, a)
    end
    pl += L2norm(vecJ, var)
    return pl
end
