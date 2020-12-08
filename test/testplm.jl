module TestPlm
using PlmIsing, Random, Test

function energy(spin,J,h)
    N =  length(h)
    enJ = 0.0
    for j in 1:N
        for i in j+1:N
            enJ -= J[i,j] * spin[i] * spin[j]
        end
    end 
    enh = 0.0
    for i in 1:N
        enh -= h[i] * spin[i]
    end
    return enJ + enh
end

function complete_dataset(N)
    Z = zeros(Int,N,2^N)
    bounds = ntuple(x->2,N)
    ctr = 0
    for i in CartesianIndices(bounds)
        ctr += 1
        for j in 1:N
            Z[j,ctr] = 2(i[j] - 1) - 1
        end
    end
    Z
end

function generateWZJh(N)
    Jasym = randn(N,N)/sqrt(N)
    h = randn(N)
    J = 0.5*(permutedims(Jasym,[2,1])+Jasym)
    for i in 1:N J[i,i] = 0.0 end
    Z = complete_dataset(N)
    W = [exp(-energy(Z[:,i],J,h)) for i in 1:2^N]
    W.= W ./sum(W)
    return W,Z,J,h
end


function test_isingplm(N;lambdaJ::Real=0.01,
    lambdaH::Real=0.01,
    epsconv::Real=1.0e-5,
    maxit::Int=1000,
    maxeval::Int=5000,
    verbose::Bool=true,
    method::Symbol=:LD_LBFGS,
    epstest::Real=1e-3)

    W,Z,J,h = generateWZJh(N)
    res=isingplm(Z,W,lambdaJ = lambdaJ,lambdaH=lambdaH,maxit=maxit,maxeval=maxeval,verbose=verbose,method=method, epsconv=epsconv)
    ΔJ =  sum(abs,(J-res.J))/length(h)^2
    Δh =  sum(abs,(h-res.H))/length(h)
    @test ΔJ < epstest
    @test Δh < epstest
    _s = 0.0
    for i in 1:2^N
        etrue = energy(Z[:,i],J,h)
        eplmi = energy(Z[:,i],res.J,res.H)
        _s += 2abs(etrue - eplmi)/abs(etrue+eplmi)
    end
    println("<Δene/ene> = ",_s/2^N, " ΔJ = $ΔJ Δh = $Δh")
end

for n in 4:12
    for it in 1:5
        test_isingplm(n;lambdaJ=0,lambdaH=0,epsconv=1e-50,verbose=false,epstest=1e-5)
    end
end

printstyled("All TestPlm passed!\n",color=:green,bold=true)
end