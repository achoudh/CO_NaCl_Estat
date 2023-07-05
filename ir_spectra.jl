# Line width
sσ(Δν) = Δν / (2.0*sqrt(2.0*log(2.0)))

# Gaussian function
gssn(ν, ν0, Δν) = exp(-(ν - ν0)^2 / (2 * sσ(Δν)^2)) / (sσ(Δν) * sqrt(2π))

# define hstatT
function hstatT(ev, eu, com0_ml, ml_in)

    h = zeros(Float64, nmols_ml, nmols_ml)

    for n1 in 1:nmols_ml

        for n2 in n1+1:nmols_ml

            rvec12 = com0_ml[n1,:] - com0_ml[n2,:] + ml_in[n1:nmols_ml:end] - ml_in[n2:nmols_ml:end] 
            rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
            rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
            rvec12 = a0_surf .* rvec12
            r12 = norm(rvec12)
            en = rvec12/r12

            force = (dot(eu[n1,:], eu[n2,:]) - 3.0*dot(en, eu[n1,:])*dot(en, eu[n2,:])) / r12^3
        
            h[n1,n1] += force
            h[n2,n2] += force
            h[n1,n2] = unit1*μ01^2*force
            h[n2,n1] = h[n1,n2]

        end

    end

    for n1 in 1:nmols_ml
        h[n1,n1] = ev[n1] + unit1*(μ11 - μ00)*μ00*h[n1,n1]
    end

    # return eigenvalues and eigenvectors
    return eigen(h)
       
end

# IR spectra
function ir_spectra(νk, x, com0_ml, Δν)

    θ = x[1+0*nmols_ml:1*nmols_ml]
    ϕ = x[1+1*nmols_ml:2*nmols_ml]
    ml_in = x[1+2*nmols_ml:5*nmols_ml]

    # Unperturbed eigenvalues
    ev = zeros(nmols_ml)
    for i in 1:nmols_ml
        ev[i] = θ[i] < 0.5*pi ? ν0[1] : ν0[2]
    end

    # Orientation of the dipole moments of each vectors
    eu = zeros(nmols_ml, 3) 
    for i in 1:nmols_ml
        eu[i,:] = [sin(θ[i]) * cos(ϕ[i]), sin(θ[i]) * sin(ϕ[i]), cos(θ[i])]
    end

    eigenvals, eigenvecs = hstatT(ev, eu, com0_ml, ml_in)
    σ = eigenvals ./ nmols_ml

    μEpda = zeros(nmols_ml)
    μEsda = zeros(nmols_ml)
    μEp = zeros(nmols_ml)
    μEs = zeros(nmols_ml)
    pl = zeros(nmols_ml,3)
    for m in 1:nmols_ml # loop over eigenvecs
        for i in 1:nmols_ml # loop over molecules
            pl[m,:] += eigenvecs[i, m]*(μ01 .* eu[i,:])
        end
        # Single domain
        μEp[m] = dot(pl[m,:],ep)^2
        μEs[m] = dot(pl[m,:],es)^2
        # Domain average
        μEpda[m] = (0.5*((pl[m,1])^2 + (pl[m,2])^2)*Tx + (pl[m,3])^2 *Tz)
        μEsda[m] = 0.5*(pl[m,1]^2 + pl[m,2]^2)*Ty
        
    end

    ipda = zeros(size(νk,1))
    isda = zeros(size(νk,1))
    ip = zeros(size(νk,1))
    is = zeros(size(νk,1))

    for (iν,ν) in enumerate(νk)
        for m in 1:nmols_ml
            gp = gssn(ν, eigenvals[m], Δν)#1.15
            gs = gssn(ν, eigenvals[m], Δν)
            ipda[iν] += unit2*σ[m]*μEpda[m] * gp
            isda[iν] += unit2*σ[m]*μEsda[m] * gs
            ip[iν] += unit2*σ[m]*μEp[m] * gp
            is[iν] += unit2*σ[m]*μEs[m] * gs
        end
    end

    return ipda, isda, ip, is
end

