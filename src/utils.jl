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

function unpackpair(J::Matrix{Float64},N::Int)
    @assert size(J,1) == 2N-1
    @assert size(J,2) == N*(N-1)>>1

    DJ = Dict{Vector{Int}, Vector{Float64}}()
    DH = [Float64[] for i=1:N]
    ctr = 0
    for i=1:N-1,j=i+1:N
        ctr += 1
        myJ = J[:,ctr]
        unpackpair!(DJ, DH, myJ, i,j,N)
    end
    Jout = zeros(N,N)
    Hout = zeros(N)
    for i=1:N-1,j=i+1:N
        Jout[i,j] = mean(DJ[[i,j]])
        Jout[j,i] = Jout[i,j]
    end
    for i=1:N
        Hout[i]=mean(DH[i])
    end
    DJ,DH,Jout,Hout
end

function unpackpair!(DJ::Dict{Vector{Int}, Vector{Float64}}, DH::Vector{Vector{Float64}},  J::Vector{Float64}, si::Int, sj::Int, N::Int)

    ctr = 0
    for i=1:si-1
        ctr += 1
        vp = sort!([i,si])
        haskey(DJ,vp) ? push!(DJ[vp],J[ctr]) : DJ[vp] = [J[ctr]]
    end
    for i=si+1:sj-1
        ctr += 1
        vp = sort!([i,si])
        haskey(DJ,vp) ? push!(DJ[vp],J[ctr]) : DJ[vp] = [J[ctr]]
    end
    for i=sj+1:N
        ctr += 1
        vp = sort!([i,si])
        haskey(DJ,vp) ? push!(DJ[vp],J[ctr]) : DJ[vp] = [J[ctr]]
    end

    for i=1:si-1
        ctr += 1
        vp = sort!([i,sj])
        haskey(DJ,vp) ? push!(DJ[vp],J[ctr]) : DJ[vp] = [J[ctr]]
    end

    for i=si+1:sj-1
        ctr += 1
        vp = sort!([i,sj])
        haskey(DJ,vp) ? push!(DJ[vp],J[ctr]) : DJ[vp] = [J[ctr]]
    end

    for i=sj+1:N
        ctr += 1
        vp = sort!([i,sj])
        haskey(DJ,vp) ? push!(DJ[vp],J[ctr]) : DJ[vp] = [J[ctr]]
    end

    ctr += 1
    vp = sort!([si,sj])
    haskey(DJ,vp) ? push!(DJ[vp],J[ctr]) : DJ[vp] = [J[ctr]]

    ctr += 1
    push!(DH[si],J[ctr])

    ctr += 1
    push!(DH[sj],J[ctr])
end

