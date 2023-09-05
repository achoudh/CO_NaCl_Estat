#############
# constants #
#############


# Two center multipole moments obtained from Meredith
mom_C::Vector{Float64} = [ 0.04135*e,  0.42424*di_au, -0.42210*qu_au]
mom_O::Vector{Float64} = [-0.04135*e, -0.26316*di_au,  0.38783*qu_au]

v1 =  0.4278e-10
v2 =  0.6649e-10
w1 = -0.8354e-10
w2 = -0.7134e-10

aij::Vector{Float64} = [613201.5, 1531479.5, 1531479.5, 5195586.0] * 1.60218e-22
bij::Vector{Float64} = [3.10, 3.49, 3.49, 4.0] * 1e10
cij::Vector{Float64} = [33096.7, 17235.6, 17235.6, 9695.3] * 1e-60 * 1.60218e-22

# k_fact::Vector{Float64} = [factorial(k) for k in 0:6]
bij_kfact::Vector{Vector{Float64}} = [(bij .^ k)/factorial(k) for k in 0:6]

#############
# functions #
#############

function ion_ion_interaction(R12::Vector{Float64}, e1z::Vector{Float64}, e2z::Vector{Float64}, 
    mom1::Vector{Float64}, mom2::Vector{Float64})

    # Determines the electrostatic interaction energy between two CO molecules

    R21 = - R12

    R = norm(R12)
    iR = 1.0 / R
    
    e12 = R12 * iR
    e21 = -e12

    # orientation vector or directional cosine

    r1z = dot(e1z, e12)
    r2z = dot(e2z, e21)
    czz = dot(e1z, e2z)

    # Multipole moments

    q1, mu1, Qu1 = mom1
    q2, mu2, Qu2 = mom2

    # interaction terms
    # 1:mu (dipole), 2:Qu(Quadrupole), 3:O(Octopole)
    # T21 means interaction between quadrupole moment of the first molecule and dipole of the second

    iR2 = iR * iR
    iR3 = iR * iR2
    iR4 = iR * iR3
    iR5 = iR * iR4
    iR6 = iR * iR5
    iR7 = iR * iR6
    
    #estat
        r1z2 = r1z * r1z
        r2z2 = r2z * r2z
        
        czz2 = czz * czz

        # charge-charge
        T00 = iR
        Vqq = q1 * q2 * T00 * eps4pi

        # charge-dipole
        T10 = iR2 * r1z
        T01 = iR2 * r2z
        Vqmu = (mu1 * q2 * T10 + q1 * mu2 * T01) * eps4pi

        #dipole-dipole
        T11 = iR3*(3*r1z*r2z + czz)
        Vmumu = mu1*mu2*T11*eps4pi

        #dipole-quadrupole
        T21 = iR4*1.5*(5.0*r1z2*r2z - r2z + 2.0*r1z*czz)
        T12 = iR4*1.5*(5.0*r2z2*r1z - r1z + 2.0*r2z*czz)
        VmuQu = (Qu1*mu2*T21 + mu1*Qu2*T12)*eps4pi

    
        #quadrupole-quadrupole
        T22=iR5*0.75*(35.0*r1z2*r2z2 - 5.0*r1z2 - 5.0*r2z2 + 20.0*r1z*r2z*czz + 2.0*czz2 + 1.0)
        VQuQu = Qu1*Qu2*T22*eps4pi

    
    V_ion_ion::Float64 = Vqq + Vqmu + Vmumu + VmuQu + VQuQu

    return V_ion_ion * joule2wn
end

function ion_ion_disp(R12::Vector{Float64}, e1z, e2z)::Float64
    o1p::Vector{Float64} = ( v2 .* e1z) 
    o2p::Vector{Float64} = ( v2 .* e2z) 
    c1p::Vector{Float64} = (-w2 .* e1z) 
    c2p::Vector{Float64} = (-w2 .* e2z) 
    
    # Distance between the atoms from two molecules
    o1o2::Float64 = norm(R12 .+ (o1p.-o2p))
    o1c2::Float64 = norm(R12 .+ (o1p.-c2p))
    c1c2::Float64 = norm(R12 .+ (c1p.-c2p))
    c1o2::Float64 = norm(R12 .+ (c1p.-o2p))

    V_disp::Float64 = 0.0

    for (i::Int64, r::Float64) in enumerate([c1c2, c1o2, o1c2, o1o2])
       fij::Float64 = 1.0
       for (k::Int64, b::Vector{Float64}) in enumerate(bij_kfact)
        fij -= b[i] * r^k * exp( - b[i]*r)
       end
       V_disp -= fij * cij[i] / r^6
    end

    return V_disp * joule2wn
end

function ion_ion_rep(R12::Vector{Float64}, e1z, e2z)::Float64
    o1p::Vector{Float64} = ( v1 .* e1z) 
    o2p::Vector{Float64} = ( v1 .* e2z) 
    c1p::Vector{Float64} = (-w1 .* e1z) 
    c2p::Vector{Float64} = (-w1 .* e2z) 
    
    # Distance between the atoms from two molecules
    o1o2::Float64 = norm(R12 .+ (o1p.-o2p))
    o1c2::Float64 = norm(R12 .+ (o1p.-c2p))
    c1c2::Float64 = norm(R12 .+ (c1p.-c2p))
    c1o2::Float64 = norm(R12 .+ (c1p.-o2p))

    V_rep::Float64 = 0.0

    for (i::Int64, r::Float64) in enumerate([c1c2, c1o2, o1c2, o1o2])
       V_rep += aij[i] * exp( - bij[i]*r)
    end

    return V_rep
end




function co_co_interaction(R12::Vector{Float64}, phi1::Float64, theta1::Float64, phi2::Float64, theta2::Float64)

 # local axis unit vectors
    stheta1, sphi1, costheta1, cosphi1 = sin(theta1), sin(phi1), cos(theta1), cos(phi1)
    stheta2, sphi2, costheta2, cosphi2 = sin(theta2), sin(phi2), cos(theta2), cos(phi2)
    
    e1z::Vector{Float64} = [stheta1 * cosphi1, stheta1 * sphi1, costheta1]
    e2z::Vector{Float64} = [stheta2 * cosphi2, stheta2 * sphi2, costheta2]

    o1p::Vector{Float64} = ( v .* e1z) 
    o2p::Vector{Float64} = ( v .* e2z) 
    c1p::Vector{Float64} = (-w .* e1z) 
    c2p::Vector{Float64} = (-w .* e2z) 
    
    # Distance between the atoms from two molecules
    o1o2::Vector{Float64} = R12 .+ (o1p.-o2p)
    o1c2::Vector{Float64} = R12 .+ (o1p.-c2p)
    c1c2::Vector{Float64} = R12 .+ (c1p.-c2p)
    c1o2::Vector{Float64} = R12 .+ (c1p.-o2p)

    V_CO_CO::Float64 = ion_ion_interaction(o1o2, e1z, e2z, mom_O, mom_O) +
                       ion_ion_interaction(o1c2, e1z, e2z, mom_O, mom_C) +
                       ion_ion_interaction(c1o2, e1z, e2z, mom_C, mom_O) +
                       ion_ion_interaction(c1c2, e1z, e2z, mom_C, mom_C) +
                       ion_ion_disp(R12, e1z, e2z) + 
                       ion_ion_rep(R12, e1z, e2z)


    return V_CO_CO
end