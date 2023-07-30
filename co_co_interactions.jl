
#########################################
# CO-CO dispersion-repulsion parameters #
#########################################

#repulsion parameters
eps_oo = 51.8037 *boltzmann # (J)
eps_cc = 31.5550 *boltzmann # (J)
# eps_co = 40.431 *boltzmann  # (J)
sig_oo = 3.0058e-10  # (m)
sig_cc = 3.56379e-10 # (m)
# sig_co = 3.2848e-10  # (m)

a_rep_oo   = 4.0*eps_oo*sig_oo^12
a_rep_d_oo = 4.0*eps_oo*sig_oo^6 / 0.78

a_rep_cc   = 4.0*eps_cc*sig_cc^12
a_rep_d_cc = 4.0*eps_cc*sig_cc^6 / 0.79

a_rep_co = sqrt( (a_rep_oo) * (a_rep_cc) )
a_rep_d_co = sqrt( (a_rep_d_oo) * (a_rep_d_cc) )


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


###############################
# Determine multipole moments #
###############################

CO_mom = multipole_R_dependence(v+w)

#############
# functions #
#############

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

function co_co_interaction(R12::Vector{Float64}, phi1::Float64, theta1::Float64, phi2::Float64, theta2::Float64)

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

    o1p = ( v .* e1z) #[ v*sin(theta1)*cos(phi1),  v*sin(theta1)*sin(phi1),  v*cos(theta1)]
    o2p = ( v .* e2z) #[ v*sin(theta2)*cos(phi2),  v*sin(theta2)*sin(phi2),  v*cos(theta2)]
    c1p = (-w .* e1z) #[-w*sin(theta1)*cos(phi1), -w*sin(theta1)*sin(phi1), -w*cos(theta1)]
    c2p = (-w .* e2z) #[-w*sin(theta2)*cos(phi2), -w*sin(theta2)*sin(phi2), -w*cos(theta2)]
    
    # Distance between the atoms from two molecules
    o1o2 = norm(R12 .+ (o1p.-o2p))
    o1c2 = norm(R12 .+ (o1p.-c2p))
    c1c2 = norm(R12 .+ (c1p.-c2p))
    o2c1 = norm(R12 .+ (c1p.-o2p))
    
    V_rep = a_rep_oo/o1o2^12 + a_rep_co/o1c2^12 + a_rep_co/o2c1^12 + a_rep_cc/c1c2^12
    
    # V_rep = a_rep_oo/o1o2^12 + a_rep_co/o1c2^12 + a_rep_co/o2c1^12 + a_rep_cc/c1c2^12
    V_disp = a_rep_d_oo/o1o2^6 + a_rep_d_co/o1c2^6 + a_rep_d_co/o2c1^6 + a_rep_d_cc/c1c2^6

    V_CO_CO::Float64 = Vqq + Vqmu + Vmumu + VmuQu + VmuO + VQuQu + VQuO + VOO

    V_CO_CO += V_rep - V_disp

    return V_CO_CO * joule2wn #Vqq, Vqmu, Vmumu, VmuQu, VmuO, VQuQu, VQuO, VOO,
end
