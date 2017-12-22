function MCgen(N::Integer, nmeasure::Integer;
               erdos::Bool = false,
               conn::Float64 = 4.0,
               nterm::Integer = 100,
               nsweepintermeasure::Integer = 10,
               p::Real=1.0,
               meanJ::Float64 = 0.9,
               h0::Float64=0.0)
    
    srand(123)

    h = h0*rand(N)
    
    if erdos
        J = er(N,conn, meanJ)
    else    
        J = InitJ(N, meanJ,p)
    end
    #   J = meanJ .* (ones(N,N) - diagm(ones(N)))./N
    s = map(x -> ((x > 0.5) ? 1.0 : -1.0), rand(N))
    H = - ( J * s  + h )

    for i = 1:nterm
        #s , H = onemcsweepslow(J,H,s,N,h)
        s , H = onemcsweep(J,H,s,N)
    end
    mag = zeros(nmeasure)
    conf = zeros(N,nmeasure)
    for i=1:nmeasure
        for j=1:nsweepintermeasure
            #s , H  = onemcsweepslow(J,H,s,N,h)
            s , H  = onemcsweep(J,H,s,N)
        end    
        mag[i] = mean(s)
        for j=1:N
            conf[j,i] = s[j]
        end
    end
   #println(maximum(abs(H + J * s  + h )))
    return mag, conf,J,h

end

function onemcsweep(J,H,s,N)
    site = randperm(N)
    @inbounds for i=1:N
        sito = site[i]
        de = -2.0 *s[sito] * H[sito]
        if de < 0
            s[sito] = -s[sito]
            scra = H[sito]
            @simd for j=1:N
                H[j] -= 2.0*J[j,sito]*s[sito]
            end
            H[sito] = scra
        elseif rand() < exp(-de)  
            s[sito] = -s[sito]
            scra = H[sito]
            @simd for j=1:N
                H[j] -= 2.0 * J[j,sito] * s[sito]
            end
            H[sito] = scra
        end
    end
    return s, H
end

function onemcsweepslow(J,H,s,N,h)
    
    site = randperm(N)
   
    for i=1:N
        sito = site[i]
        de = -2 * s[sito] * H[sito]
        if de < 0
            s[sito] = -s[sito]
            H = -(J*s + h)
        elseif rand() < exp(-de)  
            s[sito] = -s[sito]
            H = -(J*s + h)           
        end
    end
    return s, H
end

function rand_exponential(mean, N, M)
    if mean <= 0.0
        error("mean must be positive")
    end
    rndexp = -mean*log.(rand(N,M))
    return rndexp
end


function InitJ(N, meanJ, p)
    Jscra = rand_exponential(meanJ,N,N) .* map(x -> x < p, rand(N,N))
    J = triu(Jscra,1)' + triu(Jscra,1)
    ## if p < 1
    ##     J = sparse(J)
    ## end
    J = J ./ ( N * p ) 
    return J
end





function er(N,c,meanJ)
    covern = c/(2N)
    J = zeros(N,N)
    for i=1:N-1
        for j=i+1:N
            if rand()<covern
                J[i,j]= 1.
                J[j,i]= 1.
            end
        end
    end
    scale!(J,meanJ/2.)
end

## function onemcsweep(J,H,s,N,neigh)
##     site = randperm(N)
##     for i=1:N
##         sito = site[i]
##         de = -2.0 *s[sito] * H[sito]
##         if de < 0
##             s[sito] = -s[sito]
##             for j=1:N
##                 if j != sito  
##                     H[j] = H[j] - 2.0 * J[j,sito]*s[sito]
##                 end
##             end           
##         elseif rand() < exp(-de)  
##             s[sito] = -s[sito]
##             for j=1:N
##                 if j != sito  
##                     H[j] = H[j] - 2.0 * J[j,sito] * s[sito]
##                 end
##             end  
##         end
##     end
##     return s, H
## end

