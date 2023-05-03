using LinearAlgebra

#%%% Molecule Surface Stone Style


function alpha_ro_exp(params::Vector{Float64}, costheta::Float64, r::Float64)::Float64
    alpha_0::Float64, alpha_1::Float64, ro_0::Float64, ro_1::Float64, ro_2::Float64 = params
    alpha_leg::Float64 = alpha_0 + alpha_1*costheta
    ro_leg::Float64 = ro_0 + ro_1*costheta + 0.5*ro_2*(3*costheta^2 - 1)
    return exp(-alpha_leg*(r-ro_leg))
end


function rep_func(alpha, ro, r)
    return exp(-alpha*(r-ro))
end

function rel_pos_co(d, theta, phi)
    return [d*sin(theta)*cos(phi), d*sin(theta)*sin(phi), d*cos(theta)]
end

function mol_surf_rep_stone(ml_in_o::Vector{Float64}, ml_in_c::Vector{Float64}, nei::Int64)
    #costheta = cos(theta)
    op::Vector{Float64} = ml_in_o #+ [v*sin(theta)*cos(phi), v*sin(theta)*sin(phi), v*costheta]
    cp::Vector{Float64} = ml_in_c #+ [-w*sin(theta)*cos(phi), -w*sin(theta)*sin(phi), -w*costheta] 
    rco::Vector{Float64} =  op - cp
    co::Float64 = norm(rco)
    b1::Vector{Int64} = [0,1,0]
    b2::Vector{Int64} = [1,0,0]
    rep::Float64 = 0.0
    disp::Float64 = 0.0
    for i = -nei:nei
        for j = -nei:nei
            r1::Vector{Float64} = [0.0, 0.0, 0.0] + i*b1 + j*b2
            r1_2::Vector{Float64} = [0.0, 0.0, -a0_NaCl/a0_surf] + i*b1 + j*b2
            rcna::Vector{Float64} = a0_surf .* r1 - cp
            rona::Vector{Float64} = a0_surf .* r1 - op
            rccl_2::Vector{Float64} = a0_surf .* r1_2 - cp
            rocl_2::Vector{Float64} = a0_surf .* r1_2 - op
            
            
            r2::Vector{Float64} = [-0.5, -0.5, 0.0] + i*b1 + j*b2
            r2_2::Vector{Float64} = [-0.5, -0.5, -a0_NaCl/a0_surf] + i*b1 + j*b2
            rccl::Vector{Float64} = a0_surf .* r2 - cp
            rocl::Vector{Float64} = a0_surf .* r2 - op
            rcna_2::Vector{Float64} = a0_surf .* r2_2 - cp
            rona_2::Vector{Float64} = a0_surf .* r2_2 - op
            
            cna::Float64, ona::Float64, cna_2::Float64, ona_2::Float64 = norm(rcna), norm(rona), norm(rcna_2), norm(rona_2)
            cthetacna::Float64, cthetaona::Float64, cthetacna_2::Float64, cthetaona_2::Float64 = dot(rco, rcna)/(co*cna), dot(rco, rona)/(co*ona),
                                                             dot(rco, rcna_2)/(co*cna_2), dot(rco, rona_2)/(co*ona_2)

            ccl::Float64, ocl::Float64, ccl_2::Float64, ocl_2::Float64 = norm(rccl), norm(rocl), norm(rccl_2), norm(rocl_2)
            cthetaccl::Float64, cthetaocl::Float64, cthetaccl_2::Float64, cthetaocl_2::Float64 = dot(rco, rccl)/(co*ccl), dot(rco, rocl)/(co*ocl),
                                                             dot(rco, rccl_2)/(co*ccl_2), dot(rco, rocl_2)/(co*ocl_2)
            cna6::Float64, ona6::Float64, ccl6::Float64, ocl6::Float64 = [1/cna^6, 1/ona^6, 1/ccl^6, 1/ocl^6]
            cna26::Float64, ona26::Float64, ccl26::Float64, ocl26::Float64 = [1/cna_2^6, 1/ona_2^6, 1/ccl_2^6, 1/ocl_2^6]                                                    
            
            rep += (alpha_ro_exp(o_na_rep, cthetaona, ona) + alpha_ro_exp(c_na_rep, cthetacna, cna) + 
                    alpha_ro_exp(o_cl_rep, cthetaocl, ocl) + alpha_ro_exp(c_cl_rep, cthetaccl, ccl) +
                    alpha_ro_exp(o_na_rep, cthetaona_2, ona_2) + alpha_ro_exp(c_na_rep, cthetacna_2, cna_2) + 
                    alpha_ro_exp(o_cl_rep, cthetaocl_2, ocl_2) + alpha_ro_exp(c_cl_rep, cthetaccl_2, ccl_2))
                                                                                  
            disp += disp_coef[1]*cna6 + disp_coef[2]*ona6 + disp_coef[3]*ccl6 + disp_coef[4]*ocl6 +
                    disp_coef[1]*cna26 + disp_coef[2]*ona26 + disp_coef[3]*ccl26 + disp_coef[4]*ocl26
        end
    end
    rep = rep*K_stone
    return rep, disp
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

