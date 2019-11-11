function pairplmdca(filename::AbstractString;
                    lambdaJ::Real=0.005,
                    lambdaH::Real=0.01,
                    epsconv::Real=1.0e-5,
                    maxit::Int=1000,
                    maxeval::Int=5000,
                    verbose::Bool=true,
                    method::Symbol=:LD_LBFGS)

    spin = readdlm(filename)
    N,M = size(spin)
    plmalg = PlmAlg(method, verbose, epsconv, maxit, maxeval)
    plmvar = PlmVar(M, N, lambdaJ, lambdaH, spin)
    DJ, DH,outJ,outH, pslike = maximize2plmdca(plmalg,plmvar)
    PairPlmOut(sdata(pslike),DJ,DH, outJ, outH)
end

function computeX0(J0::Matrix{Float64}, H0::Vector{Float64}, si::Int, sj::Int)
    @assert si < sj
    N = length(H0)
    x0 = zeros(2N-1)
    ctr = 0
    for i=1:si-1
        ctr += 1
        x0[ctr] = J0[i,si]
    end
    for i=si+1:sj-1
        ctr += 1
        x0[ctr] = J0[i,si]
    end
    for i=sj+1:N
        ctr += 1
        x0[ctr] = J0[i,si]
    end

    for i=1:si-1
        ctr += 1
        x0[ctr] = J0[i,sj]
    end
    for i=si+1:sj-1
        ctr += 1
        x0[ctr] = J0[i,sj]
    end
    for i=sj+1:N
        ctr += 1
        x0[ctr] = J0[i,sj]
    end

    ctr += 1
    x0[ctr] = J0[si,sj]

    ctr += 1
    x0[ctr] = H0[si]

    ctr += 1
    x0[ctr] = H0[sj]

    @assert ctr == 2N-1
    return x0
end

function pairoptimfunwrapper(x::Vector, g::Vector, si,sj, var)
    g === nothing && (g = zeros(Float64, length(x)))
    return pairplmsitegrad!(x, g, si, sj, var)
end

function maximize2plmdca(alg::PlmAlg, var::PlmVar)

    @extract var M N lambdaJ lambdaH spin
    @extract alg method verbose epsconv maxit maxeval

    vecps = SharedArray{Float64}(N*(N-1)>>1)
    ctr = 0
    pairs = Array{NTuple{2,Int},1}()
    for i=1:N-1,j=i+1:N push!(pairs,(i,j)) end
    J0,H0=fill(0.0,N,N),fill(0.0,N)
    Jscra = @distributed hcat for p in pairs
        si = p[1]
        sj = p[2]
        ctr += 1
#        f(x::Vector,g::Vector)=pairplmsitegrad!(x,g,si,sj,var)
        x0 = computeX0(J0,H0,si,sj)
        opt = Opt(method,2N-1)
        ftol_abs!(opt, alg.epsconv)
        maxeval!(opt, maxit)
#        max_objective!(opt,f)
        max_objective!(opt,(x,g)->pairoptimfunwrapper(x,g,si,sj,var))
        elapsedtime = @elapsed (maxf, maxJH, ret) = optimize(opt,x0)
        alg.verbose && @printf("pair %d %d\t pll = %.4f\t time(s) = %.4f\t", si, sj, maxf, elapsedtime)
        alg.verbose && println("exit status = $ret")
        vecps[ctr] = maxf
        maxJH
    end
    DJ,DH,outJ,outH = unpackpair(Jscra,N)
    return DJ,DH,outJ,outH, vecps
end

function computeHpair(vecJ::Vector{Float64}, si::Int, sj::Int, a::Int, spin::DenseArray{Float64,2}, N::Int)
    @inbounds begin
        Hi = vecJ[end-1]
        Hj = vecJ[end]
        ctr = 0
        # begin hi
        for i=1:si-1
            ctr += 1
            Hi += vecJ[ctr] * spin[i,a]
        end
        for i=si+1:sj-1
            ctr += 1
            Hi += vecJ[ctr] * spin[i,a]
        end
        for i=sj+1:N
            ctr += 1
            Hi += vecJ[ctr] * spin[i,a]
        end

# end Hi

        for i=1:si-1
            ctr += 1
            Hj += vecJ[ctr] * spin[i,a]
        end
        for i=si+1:sj-1
            ctr += 1
            Hj += vecJ[ctr] * spin[i,a]
        end
        for i=sj+1:N
            ctr += 1
            Hj += vecJ[ctr] * spin[i,a]
        end
    end
    return Hi,Hj
end

function L2norm_pair(v::Vector, var::PlmVar)
    @extract var N lambdaJ lambdaH
    mysum = 0
    for i=1:2N-3
        mysum += v[i]^2
    end
    return lambdaJ*mysum + lambdaH*(v[end-1]^2 + v[end]^2)
end

@inline function computelogZ(Hi::Float64, Jij::Float64, Hj::Float64)
    Lcosh = log(cosh(Hi)) + log(cosh(Jij)) + log(cosh(Hj))
    Ltanh = log(1.0 + tanh(Hi)* tanh(Jij) * tanh(Hj))
    return Lcosh + Ltanh + log(4.0)
end

function gradpair!(grad::Vector{Float64}, spin::Matrix{Float64}, N::Int, Hi::Float64, Jij::Float64, Hj::Float64, si::Int, sj::Int, a::Int)
    ctr = 0

    thHi  = tanh(Hi)
    thHj  = tanh(Hj)
    thJij = tanh(Jij)



    den   = 1.0 + thHi*thJij*thHj
    num1  = thHi + thJij * thHj
    num2  = thHj + thJij * thHi
    frac1 = num1 / den
    frac2 = num2 / den


    @inbounds begin
        for i=1:si-1
            ctr += 1
            grad[ctr] += spin[si,a]*spin[i,a] - spin[i,a]*frac1
        end
        for i=si+1:sj-1
            ctr += 1
            grad[ctr] += spin[si,a]*spin[i,a] - spin[i,a]*frac1
        end
        for i=sj+1:N
            ctr += 1
        grad[ctr] += spin[si,a]*spin[i,a] - spin[i,a]*frac1
        end

        for i=1:si-1
            ctr += 1
            grad[ctr] += spin[sj,a]*spin[i,a] - spin[i,a]*frac2
        end
        for i=si+1:sj-1
            ctr += 1
            grad[ctr] += spin[sj,a]*spin[i,a] - spin[i,a]*frac2
        end
        for i=sj+1:N
            ctr += 1
            grad[ctr] += spin[sj,a]*spin[i,a] - spin[i,a]*frac2
        end

        grad[end-2] += spin[si,a]*spin[sj,a] - (thJij + thHi*thHj)/den
        grad[end-1] += spin[si,a] - frac1
        grad[end  ] += spin[sj,a] - frac2
    end
end
function pairplmsitegrad!(vecJ::Vector{Float64}, grad::Vector{Float64}, si::Int, sj::Int, var::PlmVar)
    @extract var M N lambdaJ lambdaH
    spin = sdata(var.spin)
    pl = 0.0
    Jij = vecJ[end-2]

    for i=1:2N-3
        grad[i] = 2.0 * lambdaJ * vecJ[i] * M
    end
    grad[end-1] = 2.0 * lambdaH * vecJ[end-1] * M
    grad[end  ] = 2.0 * lambdaH * vecJ[end  ] * M


    @inbounds for a=1:M
        Hi, Hj = computeHpair(vecJ, si, sj, a, spin, N)
        pl += Hi * spin[si,a] + Jij * spin[si,a] * spin[sj,a] + Hj * spin[sj,a] - computelogZ(Hi,Jij,Hj)
        gradpair!(grad,spin,N, Hi,Jij,Hj, si, sj, a)
    end

    scale!(grad,1.0/M)
    pl /= M
    pl += L2norm_pair(vecJ, var)

    return pl
end
