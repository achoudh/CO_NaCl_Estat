using LinearAlgebra

#%%% Molecule Surface Stone Style


function alpha_ro_exp(params::Vector{Float64}, costheta::Float64, r::Float64)::Float64
    alpha_0::Float64, alpha_1::Float64, ro_0::Float64, ro_1::Float64, ro_2::Float64 = params
    alpha_leg::Float64 = alpha_0 + alpha_1*costheta
    ro_leg::Float64 = ro_0 + ro_1*costheta + 0.5*ro_2*(3*costheta^2 - 1)
    return exp(-alpha_leg*(r-ro_leg))
end


function mol_surf_rep_stone(op::Vector{Float64}, cp::Vector{Float64}, nei::Int64)

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
    return repulsion*K_stone, dispersion
end


function surf_pot_stone(com::Vector{Float64}, costheta::Float64, mom::Vector{Float64})::NTuple{5, Float64}
    #costheta = cos(theta)
    
    x::Float64, y::Float64, z::Float64 = com/a0_NaCl
    # A::Float64 = eps_NaCl/a0_NaCl
    out = zeros(Float64,5)
    for i::Int64 in 1:nlm
        l::Int64, m::Int64 = lodd[i], modd[i]
        # lm::Float64 = sqrt(l^2 + m^2)
        # lma::Float64 = lm*(-2*pi/a0_NaCl)
        ret::Float64 = ret_c[i] * exp(-2*pi*z*lm[i]) * cos(2*pi*((l*x)+(m*y)-(l+m)/4))
        ret_z::Float64 = -ret*lma[i]
        ret_zz::Float64 = ret_z*lma[i]
        ret_zzz::Float64 = ret_zz*lma[i]
        ret_zzzz::Float64 = ret_zzz*lma[i]
        out += [ret, ret_z, ret_zz, ret_zzz, ret_zzzz]
    end
    
    return mom[1]*out[1], -costheta*mom[2]*out[2], -(3*(costheta)^2-1)*mom[3]*out[3]/4,
         -(5*(costheta)^3-3*costheta)*mom[4]*out[4]/12, -(35*(costheta)^4-30*costheta^2+3)*mom[5]*out[5]/192
end
    

function mol_surf_attr_stone(com_o::Vector{Float64}, com_c::Vector{Float64}, com_bc::Vector{Float64}, costheta::Float64)::Float64
    op::Vector{Float64} = rot*com_o #+ [v*sin(theta)*cos(phi), v*sin(theta)*sin(phi), v*cos(theta)]
    cp::Vector{Float64} = rot*com_c #+ [-w*sin(theta)*cos(phi), -w*sin(theta)*sin(phi), -w*cos(theta)]
    bcp::Vector{Float64} = rot*com_bc #+ [bc*sin(theta)*cos(phi), bc*sin(theta)*sin(phi), bc*cos(theta)]
    
    return sum(surf_pot_stone(op, costheta, mom_O)) + sum(surf_pot_stone(cp, costheta, mom_C)) + sum(surf_pot_stone(bcp, costheta, mom_BC))
end

