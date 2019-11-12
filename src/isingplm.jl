function isingplm(spin::Matrix{Int};
                    lambdaJ::Real=0.01,
                    lambdaH::Real=0.01,
                    epsconv::Real=1.0e-5,
                    maxit::Int=1000,
                    maxeval::Int=5000,
                    verbose::Bool=true,
                    method::Symbol=:LD_LBFGS)


    N,M = size(spin)
    plmalg = PlmAlg(method, verbose, epsconv, maxit, maxeval)
    plmvar = PlmVar(M, N, lambdaJ, lambdaH, spin)
    outJ, outh, vecps=maximizeplmdca(plmalg,plmvar)
    return PlmOut(sdata(vecps),outJ,outh)

end

function isingplm(filename::AbstractString; kwds...)
    spin = readdlm(filename,Int)
    isingplm(spin; kwds ...)
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
    ctr = 0
    # begin hi
    Hi = vecJ[end] # the i-th field
    for i=1:site-1
        ctr += 1
        Hi += vecJ[ctr] * spin[i,a]
    end
    for i=site+1:N
        ctr += 1
        Hi += vecJ[ctr] * spin[i,a]
    end
    return Hi
end

function L2norm(v::Vector, var::PlmVar)

    #pl = pl0 - lambdaJ |J|₂ - lambdaH |H|₂
    @extract var N lambdaJ lambdaH
    mysum = 0
    for i=1:N-1
        mysum += v[i]^2
    end
    return -lambdaJ*mysum - lambdaH*v[end]^2
end

function grad!(grad::Vector{Float64}, spin, N::Int, Hi::Float64, site::Int, a::Int)

    exponent = -2.0*Hi*spin[site,a]

    factor = 1.0
    if exponent < 500.0
        factor = 2.0 * spin[site,a] * exp(exponent) / (1.0 + exp(exponent))
    end

    ctr = 0
    for i=1:site-1
        ctr += 1
        grad[ctr] += spin[i,a] * factor
    end
    for i=site+1:N
        ctr += 1
        grad[ctr] += spin[i,a] * factor
    end
    grad[end] += factor
end

function plmsitegrad!(vecJ::Vector{Float64}, grad::Vector{Float64}, site::Int, var::PlmVar)
    @extract var M N lambdaJ lambdaH

    spin = sdata(var.spin)
    pl = 0.0

    for i in 1:N-1
        grad[i] = -2 * lambdaJ * vecJ[i] * M
    end
    grad[N] += -2 * lambdaH * vecJ[N] * M

    @inbounds for a=1:M
        Hi = computeH(vecJ, site, spin, a, N)
        exponent = -2.0*spin[site,a]*Hi
        plterm = -exponent
        if exponent < 700.0
            plterm = -log(1.0 +exp(-2.0*spin[site,a]*Hi))
        end

        pl += plterm
        grad!(grad, spin, N, Hi, site, a)
    end
    grad .=  grad ./ M
    pl /= M
    pl += L2norm(vecJ, var)
    return pl
end
