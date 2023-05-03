using LinearAlgebra

function monolayer(thetamono::Vector{Float64}, phimono::Vector{Float64}, z_ml::Float64)
    n = 1
    com::Matrix{Float64} = zeros(nmols_uc*nx*ny, 3)
    bondlength::Vector{Float64} = zeros(nmols_uc*nx*ny)
    theta::Vector{Float64} = zeros(nmols_uc*nx*ny)
    phi::Vector{Float64} = zeros(nmols_uc*nx*ny)

    xC_unitcell = zeros(Float64, nmols_uc, 3)
    #xO_unitcell = zeros(nmols_uc, 3)
    xC_unitcell = ucell_ml(0.0, thetamono, phimono)

    for i in 0:nx-1
        for j in 0:ny-1
            com[n:n+nmols_uc-1,1] .= xC_unitcell[:,1] .+ 2*i
            com[n:n+nmols_uc-1,2] .= xC_unitcell[:,2] .+ 2*j
            com[n:n+nmols_uc-1,3] .= xC_unitcell[:,3] .+ z_ml
            n += nmols_uc
        end
    end

    for i in 1:nmols_uc
        bondlength[i] = (v + w) #/a0_surf
        theta[i] = thetamono[i]
        phi[i] = phimono[i]
    end

    i = nx*ny
    j = nmols_uc*i
    bondlength = vec(repeat(bondlength[1:nmols_uc], outer=(1,i)))
    theta = vec(repeat(theta[1:nmols_uc], outer=(1,i)))
    phi = vec(repeat(phi[1:nmols_uc], outer=(1,i)))
    #phi = reshape(repeat(phi[1:nmols_uc], outer=(1,i)), j)

    return com, bondlength, phi, theta
end

function overlayer(thetaol::Vector{Float64}, phiol::Vector{Float64}, z_ml::Float64, dz::Float64)
    n = 1
    com::Matrix{Float64} = zeros(nmols_ucol*nx*ny*nz, 3)
    bondlength::Vector{Float64} = zeros(nmols_ucol*nx*ny*nz)
    theta::Vector{Float64} = zeros(nmols_ucol*nx*ny*nz)
    phi::Vector{Float64} = zeros(nmols_ucol*nx*ny*nz)

    xC_unitcell = zeros(Float64, nmols_ucol,3)
    xO_unitcell = zeros(Float64, nmols_ucol,3)
    xC_unitcell = ucell_alpha(0.0, thetaol, phiol)

    for i in 0:nx-1
        for j in 0:ny-1
            for k in 0:nz-1
                com[n:n+nmols_ucol-1,1] .= xC_unitcell[:,1] .+ 2*i
                com[n:n+nmols_ucol-1,2] .= xC_unitcell[:,2] .+ 2*j
                com[n:n+nmols_ucol-1,3] .= xC_unitcell[:,3] .+ (a0_CO/a0_surf) .* (k + dz) .+ z_ml
                n += nmols_ucol
            end
        end
    end

    for i in 1:nmols_ucol
        bondlength[i] = (v + w) #/a0_surf
        theta[i] = thetaol[i]
        phi[i] = phiol[i]
    end

    i = nx*ny*nz
    j = nmols_ucol*i
    bondlength = vec(repeat(bondlength[1:nmols_ucol], outer=(1,i)))#, j)
    theta = vec(repeat(theta[1:nmols_ucol], outer=(1,i)))#, j)
    phi = vec(repeat(phi[1:nmols_ucol], outer=(1,i)))#, j)

    return com, bondlength, phi, theta
end

function ucell_alpha(v::Float64, thetaol::Vector{Float64}, phiol::Vector{Float64})
    n_mols_ucol::Int64 = length(thetaol)
    coords = zeros(Float64,n_mols_ucol,3)
    coxyz = zeros(Float64,n_mols_ucol,3)

    coxyz[:,1] = v .* sin.(thetaol) .* cos.(phiol)
    coxyz[:,2] = v .* sin.(thetaol) .* sin.(phiol)
    coxyz[:,3] = v .* cos.(thetaol)
    
    coords[1,:] = [0.5 + coxyz[1,1], 0.5 + coxyz[1,2], coxyz[1,3]]
    coords[2,:] = [1.5 + coxyz[2,1], 0.5 + coxyz[2,2], coxyz[2,3]]
    coords[3,:] = [1.5 + coxyz[3,1], 1.5 + coxyz[3,2], coxyz[3,3]]
    coords[4,:] = [0.5 + coxyz[4,1], 1.5 + coxyz[4,2], coxyz[4,3]]
    coords[5,:] = [0.0 + coxyz[5,1], 0.0 + coxyz[5,2], (a0_CO/a0_surf) * 0.5 + coxyz[5,3]]
    coords[6,:] = [1.0 + coxyz[6,1], 0.0 + coxyz[6,2], (a0_CO/a0_surf) * 0.5 + coxyz[6,3]]
    coords[7,:] = [1.0 + coxyz[7,1], 1.0 + coxyz[7,2], (a0_CO/a0_surf) * 0.5 + coxyz[7,3]]
    coords[8,:] = [0.0 + coxyz[8,1], 1.0 + coxyz[8,2], (a0_CO/a0_surf) * 0.5 + coxyz[8,3]]
    return coords
end

function ucell_ml(v::Float64, thetamono::Vector{Float64}, phimono::Vector{Float64})
    coxyz = zeros(Float64,length(phimono),3)
    coords = zeros(Float64,length(phimono),3)

    coxyz[:,1] .= v .* sin.(thetamono) .* cos.(phimono)
    coxyz[:,2] .= v .* sin.(thetamono) .* sin.(phimono)
    coxyz[:,3] .= v .* cos.(thetamono)
    coords[1,:] .= [0.0 + coxyz[1,1], 0.0 + coxyz[1,2], coxyz[1,3]]
    coords[2,:] .= [1.0 + coxyz[2,1], 0.0 + coxyz[2,2], coxyz[2,3]]
    coords[3,:] .= [1.0 + coxyz[3,1], 1.0 + coxyz[3,2], coxyz[3,3]]
    coords[4,:] .= [0.0 + coxyz[4,1], 1.0 + coxyz[4,2], coxyz[4,3]]

    return coords
end

function ucell_ml_2x2(v::Float64, coords::Matrix{Float64})
    coords[:,1] .= [v, 0.5 - v, 0.5 + v]
    coords[:,2] .= [0.5 - v, -v, 0.5 + v]
end

function kronecker_delta(a::Int, b::Int)::Int
    if a == b
        return 1
    else
        return 0
    end
end


       
