function unpack(J::Matrix{Float64})

    N,N1 = size(J)
    N != N1 && error("not a square matrix")

    Jret = fill(0.0, N, N)
    hret = fill(0.0, N)

    
    for i in 1:N
        hret[i] = J[N,i]
        ctr = 0
        for j in 1:i-1
            ctr += 1
            Jret[j,i] = J[ctr,i]
        end
        for j in i+1:N
            ctr += 1
            Jret[j,i] = J[ctr,i]
        end
    end
    Jret, hret
end
