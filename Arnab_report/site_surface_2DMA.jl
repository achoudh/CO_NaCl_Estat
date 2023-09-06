using LinearAlgebra

######################################
# CO-NaCl electrostatic              #
# 2 centered DMA from Hoang et al.   #
# Buckingham style field interaction #
######################################


#############
# constants #
#############

# Two center multipole moments obtained from Meredith
mom_C::Vector{Float64} = [ 0.04135*e,  0.42424*di_au, -0.42210*qu_au]
mom_O::Vector{Float64} = [-0.04135*e, -0.26316*di_au,  0.38783*qu_au]


# All constants below are reported as alpha0, alpha1, ro0, ro1, ro2
c_na_rep::Vector{Float64} = [4.106e10, 2.456e-10, -0.795e-10, -0.239e-10]
c_cl_rep::Vector{Float64} = [3.260e10, 2.989e-10, -0.953e-10, -0.267e-10]
o_na_rep::Vector{Float64} = [5.085e10, 2.942e-10, -0.374e-10,  0.101e-10]
o_cl_rep::Vector{Float64} = [3.880e10, 3.616e-10, -0.387e-10,  0.129e-10]
rep_coeffs = [ [c_na_rep, o_na_rep], [c_cl_rep, o_cl_rep] ]
K_stone::Float64 = 27.25 * 1.60218e-22

# Dispersion coefficients
# [ [C-Na, O-Na], [C-Cl, O-Cl] ]
disp_coef::Matrix{Float64} = [ [3774 2889]; [38740 29655] ]*1e-60 * 1.60218e-22 ## meV to J


#############
# functions #
#############


function multipole_components(e1z::Array{Float64, 1}, mom::Vector{Float64})
    # Calculate components of the multipole moments
    # mu, Qu, o -> are the scalar values of dipole moment, quadrapole moment and octapole moment respectively
    # All the terms in Cartesian Coordinate are calculated from (Buckingham, 1959)
    
    charge::Float64 = 0.0
    dipole::Float64 = 0.0
    quadrupole::Float64 = 0.0

    Q0::Float64 = 0.0
    Q1::Array{Float64, 1} = zeros(3)
    Q2::Array{Float64, 2} = zeros(3, 3)
    
    # Call multipole_R_dependence subroutine to get multipole moments
    charge, dipole, quadrupole = mom

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

    return Q0, Q1, Q2
end

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

function site_surface_interaction(com::Vector{Float64}, unit_vec::Vector{Float64}, mom::Vector{Float64})
    
    # Return variables
    V_CO_NaCl = 0.0
    V_Q0_NaCl = 0.0
    V_Q1_NaCl = zeros(3)
    V_Q2_NaCl = zeros(3, 3)
    
    com1::Vector{Float64} = rot*com
    
    # Orientation unit vector
    # stheta = sin(theta1)
    e1z = unit_vec #[stheta*cos(phi1), stheta*sin(phi1), cos(theta1)]

    #println(e1z)
    # Calculate multipole moments
    Q0, Q1, Q2 = multipole_components(e1z, mom)

    # Calculate electric field and its derivatives
    potfactor = eps_NaCl/a0_NaCl
    phi, phia, phiab = zeros(1), zeros(3), zeros(3,3)

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

 

    # Total CO-NaCl interaction energy is
    V_CO_NaCl = V_Q0_NaCl + V_Q1_NaCl + V_Q2_NaCl

    return V_CO_NaCl, V_Q0_NaCl, V_Q1_NaCl, V_Q2_NaCl
    
end

function mol_surf_attr_2DMA_tensor(com_o::Vector{Float64}, com_c::Vector{Float64}, unit_vec::Vector{Float64})::Float64
    
    return  site_surface_interaction(com_o, unit_vec, mom_O)[1] + 
            site_surface_interaction(com_c, unit_vec, mom_C)[1] 
end


####################################
# CO-NaCl Repulsion and Dispersion #
####################################

function alpha_ro_exp(params::Vector{Float64}, costheta::Float64, r::Float64)::Float64
    alpha_0::Float64, ro_0::Float64, ro_1::Float64, ro_2::Float64 = params
    ro_leg::Float64 = ro_0 + ro_1*costheta + 0.5*ro_2*(3*costheta^2 - 1)
    return exp(-alpha_0*(r-ro_leg))
end

function mol_surf_rep_2DMA(op::Vector{Float64}, cp::Vector{Float64}, nei::Int64)

    mco::Vector{Vector{Float64}} = [cp, op]
    rco::Vector{Float64} =  op - cp
    dco::Float64 = norm(rco)
    uco::Vector{Float64} = rco/dco
    repulsion::Float64 = 0.0
    dispersion::Float64 = 0.0
    r_nei = zeros(Float64, 3)
    Δr = zeros(Float64, 3)
    dist ::Float64 = 0.0
    cosθ ::Float64 = 0.0
    #i :: Int16,j :: Int16,layer :: Int16,ion :: Int16,atom :: Int16 = 1,1,1,1,1

    for i:: Int16 = -nei:nei
    for j:: Int16 = -nei:nei
        r_nei = i*b1 + j*b2

       for layer:: Int16 in 1:nl_surf
       for ion:: Int16 in 1:2
       for atom:: Int16 in 1:2

            Δr =  pos_surf[layer][:,ion] + r_nei - mco[atom]
            dist = norm(Δr)
            cosθ = uco⋅Δr/dist

            repulsion += alpha_ro_exp(rep_coeffs[ion][atom], cosθ, dist)
            dispersion += disp_coef[ion,atom]/dist^6

       end
       end
       end

    end
    end

    return repulsion*K_stone - dispersion
end

