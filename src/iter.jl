function isingplmdca(filename::AbstractString;
                    lambdaJ::Real=0.01,
                    lambdaH::Real=0.01,
                    epsconv::Real=1.0e-5,
                    maxit::Int=1000,
                    maxeval::Int=5000,
                    verbose::Bool=true,
                    method::Symbol=:LD_LBFGS)


    spin = readdlm(filename,Int)

    N,M = size(spin)
    plmalg = PlmAlg(method, verbose, epsconv, maxit, maxeval)
    plmvar = PlmVar(M, N, lambdaJ, lambdaH, spin)
    maximizeplmdca(plmalg,plmvar)
    
#    DJ, DH,outJ,outH, pslike = maximizeplmdca(plmalg,plmvar, J0, H0)    
#    PlmOut(sdata(pslike),DJ,DH, outJ, outH)
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
    Jscra = @parallel hcat for i in 1:N
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
    Hi = vecJ[end] # the i-th field
    ctr = 0
    # begin hi
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
    @extract var N lambdaJ lambdaH
    mysum = 0
    for i=1:N-1
        mysum += v[i]^2
    end       
    return lambdaJ*mysum + lambdaH*v[end]^2 
end


function grad!(grad::Vector{Float64}, spin, N::Int, Hi::Float64, site::Int, a::Int)

    factor = 2.0 * spin[site,a] * exp(-2.0*Hi*spin[site,a]) / (1.0 + exp(-2Hi*spin[site,a]))    
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
    @inbounds for a=1:M
        Hi = computeH(vecJ, site, spin, a, N)        
        pl += -log(1.0 +exp(-2.0*spin[site,a]*Hi))
        grad!(grad, spin, N, Hi, site, a)
        for i in 1:N-1 grad[i] += 2 * lambdaJ * vecJ[i] end
        grad[N] += 2 * lambdaJ * vecJ[N] 
    end
    scale!(grad, 1.0/M )
    pl /= M
    pl += L2norm(vecJ, var)
    return pl
end
