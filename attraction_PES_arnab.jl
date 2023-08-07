#######################################
# CO R-dependent multipole components #
#######################################

# dipole
a10::Float64 =  3.63272
a11::Float64 = -7.28409
a12::Float64 =  7.5581
a13::Float64 = -4.28056
a14::Float64 =  1.30636
a15::Float64 = -0.215764
a16::Float64 =  0.0182316
a17::Float64 = -0.000619188
# quadrupole
a20::Float64 = -14.5407
a21::Float64 =  28.2314
a22::Float64 = -25.6766
a23::Float64 =  11.9892
a24::Float64 =  -2.99969
a25::Float64 =   0.404491
a26::Float64 =  -0.0273732
a27::Float64 =   0.00071271
# octapole
a30::Float64 =  6.45176
a31::Float64 = -10.3863
a32::Float64 =  2.47699
a33::Float64 =  5.54757
a34::Float64 = -3.83351
a35::Float64 =  0.969093
a36::Float64 = -0.108525
a37::Float64 =  0.00450608

function multipole_R_dependence(R::Float64)#, charge::Float64, dipole::Float64, quadrupole::Float64, octapole::Float64)
    
    # Calculates multipole moments of the molecules
    # R : Intermolecular distance
    
    # All the parameters are calculated in the Mathematica Notebook "P:\Research\CO-NaCl\Theory\Electrostatic model\MM_R_fit\R_MM_fit.nb"
    # Data is collected from Jascha.

    R1 = R*m2au
    R1_2 = R1*R1
    R1_3 = R1*R1_2
    R1_4 = R1*R1_3
    R1_5 = R1*R1_4
    R1_6 = R1*R1_5
    R1_7 = R1*R1_6

    charge = 0.0
    dipole     = a10 + a11*R1 + a12*R1_2 + a13*R1_3 + a14*R1_4 + a15*R1_5 + a16*R1_6 + a17*R1_7
    quadrupole = a20 + a21*R1 + a22*R1_2 + a23*R1_3 + a24*R1_4 + a25*R1_5 + a26*R1_6 + a27*R1_7
    octapole   = a30 + a31*R1 + a32*R1_2 + a33*R1_3 + a34*R1_4 + a35*R1_5 + a36*R1_6 + a37*R1_7

    # Conversion to SI
    dipole     = di_au*dipole
    quadrupole = qu_au*quadrupole
    octapole   = oc_au*octapole
    
    return charge, dipole, quadrupole, octapole
end

###############################
# Determine multipole moments #
###############################

CO_mom = multipole_R_dependence(v+w)

#########################################
# CO-NaCl Attraction: multipoles at COM #
#########################################

function surf_pot_arnab(com::Vector{Float64}, costheta::Float64, mom::Vector{Float64})::NTuple{4, Float64}

    x, y, z = com/a0_NaCl
    ret = zeros(Float64, 5)
    for i::Int64 in 1:nlm
        l::Int64, m::Int64 = lodd[i], modd[i]
        temp = ret_c[i] * exp(-2*pi*z*lm[i]) * cos(2*pi*((l*x)+(m*y)-(l+m)/4))
        ret[1] += temp
        for j in 2:4  # upto 4th moment
           temp *= lma[i]
           ret[j] -= temp 
        end
    end
    
    return mom[1]*ret[1], -costheta*mom[2]*ret[2], -(3*(costheta)^2-1)*mom[3]*ret[3]/4,
         -(5*(costheta)^3-3*costheta)*mom[4]*ret[4]/12
end
    

function mol_surf_attr_arnab(com::Vector{Float64}, costheta::Float64)::Float64
    comp::Vector{Float64} = rot*com
    mom = collect(multipole_R_dependence(v+w))
    return sum(surf_pot_arnab(comp, costheta, mom))
end

##########################################
# CO-CO electrostatic: multipoles at COM #
##########################################

function co_co_interaction(R12::Vector{Float64}, phi1::Float64, theta1::Float64, phi2::Float64, theta2::Float64,
    rep_on::Int64=1, disp_on::Int64=1)

    # Determines the electrostatic interaction energy between two CO molecules
    # using electrostatic interaction, dispersion and repulsion

    R21 = - R12

    R = norm(R12)
    iR = 1.0 / R

    e12 = R12 * iR
    e21 = -e12

    # local axis unit vectors
    stheta1, sphi1, costheta1, cosphi1 = sin(theta1), sin(phi1), cos(theta1), cos(phi1)
    stheta2, sphi2, costheta2, cosphi2 = sin(theta2), sin(phi2), cos(theta2), cos(phi2)

    e1z = [stheta1 * cosphi1, stheta1 * sphi1, costheta1]
    e2z = [stheta2 * cosphi2, stheta2 * sphi2, costheta2]

    # orientation vector or directional cosine

    r1z = dot(e1z, e12)
    r2z = dot(e2z, e21)
    czz = dot(e1z, e2z)

    # Multipole moments

    q1, mu1, Qu1, O1 = CO_mom
    q2, mu2, Qu2, O2 = CO_mom

    # interaction terms
    # 1:mu (dipole), 2:Qu(Quadrupole), 3:O(Octopole)
    # T21 means interaction between quadrupole moment of the first molecule and dipole of the second

    iR2 = iR * iR
    iR3 = iR * iR2
    iR4 = iR * iR3
    iR5 = iR * iR4
    iR6 = iR * iR5
    iR7 = iR * iR6

    r1z2 = r1z * r1z
    r2z2 = r2z * r2z
    r1z3 = r1z * r1z2
    r2z3 = r2z * r2z2
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

    #dipole-octopole
    T31 = iR5*0.5*( 5.0*(7.0*r1z3 - 3.0*r1z)*r2z + 3.0*(5.0*r1z2 - 1.0)*czz )
    T13 = iR5*0.5*( 5.0*(7.0*r2z3 - 3.0*r2z)*r1z + 3.0*(5.0*r2z2 - 1.0)*czz )
    VmuO = (O1*mu2*T31 + mu1*O2*T13)*eps4pi

    #quadrupole-quadrupole
    T22=iR5*0.75*(35.0*r1z2*r2z2 - 5.0*r1z2 - 5.0*r2z2 + 20.0*r1z*r2z*czz + 2.0*czz2 + 1.0)
    VQuQu = Qu1*Qu2*T22*eps4pi

    #quadrupole-octopole
    T32 = iR6*1.25*(21.0*(3.0*r1z3 - r1z)*r2z2 + 6.0*(7.0*r1z2 - 1.0)*r2z*czz - 7.0*r1z3 
    + 3.0*r1z+6.0*r1z*czz2)
    T23 = iR6*1.25*(21.0*(3.0*r2z3 - r2z)*r1z2 + 6.0*(7.0*r2z2 - 1.0)*r1z*czz - 7.0*r2z3
    + 3.0*r2z+6.0*r2z*czz2)
    VQuO = (O1*Qu2*T32 + Qu1*O2*T23)*eps4pi

    #octopole-octopole
    T33 = iR7*1.25*(21.0*(11.0*r1z3 - 3.0*r1z)*r2z3 - 21.0*(3.0*r1z3 - r1z)*r2z
    + (189.0*r1z2*r2z2 - 21.0*r1z2 - 21.0*r2z2 + 3.0)*czz + 42.0*r1z*r2z*czz2 + 2.0*czz2*czz)
    VOO = O1*O2*T33*eps4pi

   
    V_CO_CO::Float64 = Vqq + Vqmu + Vmumu + VmuQu + VmuO + VQuQu + VQuO + VOO

    return V_CO_CO * joule2wn #Vqq, Vqmu, Vmumu, VmuQu, VmuO, VQuQu, VQuO, VOO,
end
