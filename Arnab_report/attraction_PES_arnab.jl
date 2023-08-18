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

function multipole_components(R::Float64, e1z::Array{Float64, 1})
    # Calculate components of the multipole moments
    # mu, Qu, o -> are the scalar values of dipole moment, quadrapole moment and octapole moment respectively
    # All the terms in Cartesian Coordinate are calculated from (Buckingham, 1959)
    
    charge::Float64 = 0.0
    dipole::Float64 = 0.0
    quadrupole::Float64 = 0.0
    octapole::Float64 = 0.0
    Q0::Float64 = 0.0
    Q1::Array{Float64, 1} = zeros(3)
    Q2::Array{Float64, 2} = zeros(3, 3)
    Q3::Array{Float64, 3} = zeros(3, 3, 3)
    
    # Call multipole_R_dependence subroutine to get multipole moments
    charge, dipole, quadrupole, octapole = multipole_R_dependence(R)

    # Calculating Cartesian components
    # charge
    Q0 = charge

    # Dipole moments
    Q1 = dipole * e1z

    # Quadrupole moments
    for i in 1:3
        for j in 1:3
            Q2[i, j] = quadrupole * (3.0 * e1z[i] * e1z[j] - (i == j)) * 0.5
        end
    end

    # Octapole moments
    for i in 1:3
        for j in 1:3
            for k in 1:3
                Q3[i, j, k] = octapole * (5.0 * e1z[i] * e1z[j] * e1z[k] - e1z[i] * (j == k) -
                                          e1z[j] * (k == i) - e1z[k] * (i == j)) * 0.5
            end
        end
    end

    return Q0, Q1, Q2, Q3
end


###############################
# Determine multipole moments #
###############################

CO_mom = multipole_R_dependence(v+w)

#########################################
# CO-NaCl Attraction: multipoles at COM #
#        Scalar multiplication          #
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

#########################################
# CO-NaCl Attraction: multipoles at COM #
#        Vector multiplication          #
#########################################
function pot_deriv_lm(com1::Vector{Float64}, dxyz::Vector{Float64}, l::Int64, m::Int64)::Float64
    # This function calculates l and m dependent terms for the partial derivatives of electric potential
    la = l/a0_NaCl
    ma = m/a0_NaCl
    sqrt_lm = sqrt(l*l + m*m)
    sqrt_lma = sqrt_lm/a0_NaCl

    if mod(dxyz[1]+dxyz[2], 4) == 0 # derivatives of the cos term with respect to x and y
        cosderiv = cos(2.0*pi*(la*com1[1] + ma*com1[2] - 0.25*(l+m))) * (2.0*pi*la)^dxyz[1] * (2.0*pi*ma)^dxyz[2]
    elseif mod(dxyz[1]+dxyz[2], 4) == 1
        cosderiv = -sin(2.0*pi*(la*com1[1] + ma*com1[2] - 0.25*(l+m))) * (2.0*pi*la)^dxyz[1] * (2.0*pi*ma)^dxyz[2]
    elseif mod(dxyz[1]+dxyz[2], 4) == 2
        cosderiv = -cos(2.0*pi*(la*com1[1] + ma*com1[2] - 0.25*(l+m))) * (2.0*pi*la)^dxyz[1] * (2.0*pi*ma)^dxyz[2]
    elseif mod(dxyz[1]+dxyz[2], 4) == 3
        cosderiv = sin(2.0*pi*(la*com1[1] + ma*com1[2] - 0.25*(l+m))) * (2.0*pi*la)^dxyz[1] * (2.0*pi*ma)^dxyz[2]
    end

    # derivatives of the exponential term with respect to z
    expderiv = (-2.0*pi*sqrt_lma)^dxyz[3] * exp(-2.0*π*com1[3]*sqrt_lma)

    # summand for a single (l,m) pair
    return (-1.0)^((l+m)/2) / sqrt_lm * cosderiv * expderiv / (1.0 + exp(-π*sqrt_lm))
end

function site_surface_interaction(com::Vector{Float64}, theta1::Float64, phi1::Float64)
    
    # Return variables
    V_CO_NaCl = 0.0
    V_Q0_NaCl = 0.0
    V_Q1_NaCl = zeros(3)
    V_Q2_NaCl = zeros(3, 3)
    V_Q3_NaCl = zeros(3, 3, 3)
    V_rep_surf = 0.0
    
    com1::Vector{Float64} = rot*com
    
    # Orientation unit vector
    stheta = sin(theta1)
    e1z = [stheta*cos(phi1), stheta*sin(phi1), cos(theta1)]
    #println(e1z)
    # Calculate multipole moments
    Q0, Q1, Q2, Q3 = multipole_components(v+w, e1z)
    #println(Q0,Q1,Q2,Q3)
    # Calculate electric field and its derivatives
    potfactor = eps_NaCl/a0_NaCl
    phi, phia, phiab, phiabc = zeros(1), zeros(3), zeros(3,3), zeros(3,3,3)
    for i = 1:nlm
        # Electric potential
        dxyz = zeros(3)
        phi[1] += potfactor*pot_deriv_lm(com1, dxyz, lodd[i], modd[i])

        # Electric field
        for k1 = 1:3
        dxyz = zeros(3)
        dxyz[k1] += 1
        phia[k1] -= potfactor*pot_deriv_lm(com1, dxyz, lodd[i], modd[i])
        end

        # Electric field gradient
        for k1 = 1:3, k2 = 1:3
        dxyz = zeros(3)
        dxyz[k1] += 1
        dxyz[k2] += 1
        phiab[k1, k2] -= potfactor*pot_deriv_lm(com1, dxyz, lodd[i], modd[i])
        end

        # Electric field hypergradient
        for k1 = 1:3, k2 = 1:3, k3 = 1:3
        dxyz = zeros(3)
        dxyz[k1] += 1
        dxyz[k2] += 1
        dxyz[k3] += 1
        phiabc[k1, k2, k3] -= potfactor*pot_deriv_lm(com1, dxyz, lodd[i], modd[i])
        end
    end

    V_Q0_NaCl = Q0*phi[1]

    # dipole-electric field interaction
    V_Q1_NaCl = 0.0
    for k1 in 1:3
        V_Q1_NaCl -= Q1[k1]*phia[k1]
    end

    # quadrupole-electric field gradient interaction
    V_Q2_NaCl = 0.0
    for k1 in 1:3, k2 in 1:3
        V_Q2_NaCl -= Q2[k1,k2]*phiab[k1,k2]/3.0
    end

    # octapole-second gradient interaction
    V_Q3_NaCl = 0.0
    for k1 in 1:3, k2 in 1:3, k3 in 1:3
        V_Q3_NaCl -= Q3[k1,k2,k3]*phiabc[k1,k2,k3]/15.0
    end

    # Total CO-NaCl interaction energy is
    V_CO_NaCl = V_Q0_NaCl + V_Q1_NaCl + V_Q2_NaCl + V_Q3_NaCl
    # println(V_CO_NaCl*joule2wn,V_Q1_NaCl*joule2wn,  V_Q2_NaCl*joule2wn, V_Q3_NaCl*joule2wn)

    return V_CO_NaCl, V_Q0_NaCl, V_Q1_NaCl, V_Q2_NaCl, V_Q3_NaCl, V_rep_surf 
    
end

##########################################
# CO-CO electrostatic: multipoles at COM #
##########################################

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

   
    V_CO_CO::Float64 = Vqq + Vqmu + Vmumu + VmuQu + VmuO + VQuQu + VQuO + VOO

    return V_CO_CO * joule2wn #Vqq, Vqmu, Vmumu, VmuQu, VmuO, VQuQu, VQuO, VOO,
end
