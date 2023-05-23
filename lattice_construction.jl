# Defines c.-of-m. and orientation of the molecules in monolayer
function monolayer(θ_uc::Vector{Float64}, ϕ_uc::Vector{Float64}, z::Float64)

    com::Matrix{Float64} = zeros( nmols_uc*nx*ny, 3)

    com_uc::Matrix{Float64} = [ [0.0 0.0]; 
                                [1.0 0.0];
                                [1.0 1.0];
                                [0.0 1.0] ]
    n = 1
    for i in 0:nx-1
        for j in 0:ny-1
            com[n:n+nmols_uc-1,1] = com_uc[:,1] .+ 2*i
            com[n:n+nmols_uc-1,2] = com_uc[:,2] .+ 2*j
            com[n:n+nmols_uc-1,3] .= z
            n += nmols_uc
        end
    end

    phi = repeat(ϕ_uc, outer=nx*ny)
    theta = repeat(θ_uc, outer=nx*ny)

    return com, phi, theta
end

# Defines c.-of-m. and orientation of the molecules in overlayer
function overlayer(θ_uc::Vector{Float64}, ϕ_uc::Vector{Float64}, z::Float64, dz::Float64)

    com = zeros(nmols_ucol*nx*ny*nz, 3)

    a = 0.5*a0_CO/a0_surf
    com_uc::Matrix{Float64} = [ [0.5 0.5 0.0]; 
                                [1.5 0.5 0.0];
                                [1.5 1.5 0.0];
                                [0.5 1.5 0.0];
                                [0.0 0.0   a];
                                [1.0 0.0   a];
                                [1.0 1.0   a];
                                [0.0 1.0   a] ]

    n = 1
    for i in 0:nx-1
        for j in 0:ny-1
            for k in 0:nz-1
                com[n:n+nmols_ucol-1,1] .= com_uc[:,1] .+ 2*i
                com[n:n+nmols_ucol-1,2] .= com_uc[:,2] .+ 2*j
                com[n:n+nmols_ucol-1,3] .= com_uc[:,3] .+ (2*a*(k + dz) + z)
                n += nmols_ucol
            end
        end
    end

    phi = repeat(ϕ_uc, outer=nx*ny*nz)
    theta = repeat(θ_uc, outer=nx*ny*nz)

    return com, phi, theta
end