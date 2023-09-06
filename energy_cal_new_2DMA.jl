##############################################
# The energy calculations are being modified #
# after finding that 3 center DMA also       #
# required tensor multiplication             #
##############################################

include("Arnab_report/site_surface_2DMA.jl")
# include("Arnab_report/hoang_interactions.jl")
include("co_co_Guo.jl")


# Total energy (intra-monolayer + intra-overlayer + monolayer-Surface + overlayer-Surface + monolayer-overlayer )
function energy(x, lattice_ml, lattice_ol, ϕ_ol, θ_ol) 

    theta_ml = x[1+0*nmols_ml:1*nmols_ml]
    phi_ml =   x[1+1*nmols_ml:2*nmols_ml]
    δr_ml =    x[1+2*nmols_ml:5*nmols_ml]
    δz_ol =    x[ndofs_ml]
    theta_ol = x[1+ndofs_ml+0*nmols_ol2:ndofs_ml+1*nmols_ol2]
    phi_ol   = x[1+ndofs_ml+1*nmols_ol2:ndofs_ml+2*nmols_ol2]
    xy_ol    = x[1+ndofs_ml+2*nmols_ol2:ndofs_ml+4*nmols_ol2]
    
    # intra-monolayer interaction

    pot_mlml = Float64(0.0)
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        rvec12 = lattice_ml[i,:] - lattice_ml[j,:] + δr_ml[i:nmols_ml:end] - δr_ml[j:nmols_ml:end] 
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlml += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end

    # intra-overlayer interaction
    # println("mol. 1","\t", "mol. 2","\t", "Distance/Å","\t \t", "Energy/cm-1")
    pot_olol = Float64(0.0)
    for i in 1:nmols_ol2-1, j in i+1:nmols_ol2
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] +[ xy_ol[i:nmols_ol2:end] - xy_ol[j:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_olol += co_co_interaction(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j])
        # println(i,"\t", j,"\t", norm(rvec12)/1e-10,"\t \t", co_co_interaction(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j]))
    end
    # println(pot_olol)
    for i in 1:nmols_ol2, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] + [xy_ol[i:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        # println(rvec12)
        pot_olol += co_co_interaction(rvec12, phi_ol[i], theta_ol[i], ϕ_ol[j], θ_ol[j])
    end

    # Monolayer-Surface interaction

    pot_mlsurf = Float64(0.0)
    for i in 1:nmols_ml

        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = [δr_ml[i], δr_ml[i+nmols_ml], δr_ml[i+2*nmols_ml] + lattice_ml[i,3]] * a0_surf
        unit_vec = [stheta*cosphi, stheta*sphi, costheta]

        ml_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        pot_mlsurf += mol_surf_attr_2DMA_tensor(ml_o, ml_c, unit_vec) +
                      mol_surf_rep_2DMA(ml_o, ml_c, 4)
    end

    # Overlayer-Surface interaction

    pot_olsurf = Float64(0.0)
    for i in 1:nmols_ol2

        rvec = lattice_ol[i,:]
        rvec[1] -= round(rvec[1])
        rvec[2] -= round(rvec[2])

        stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
        rvec += [xy_ol[i], xy_ol[i+nmols_ol2], δz_ol] 
        rvec *= a0_surf
        unit_vec = [stheta*cosphi, stheta*sphi, costheta]

        ol_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ol_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ol_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        pot_olsurf += mol_surf_attr_2DMA_tensor(ol_o, ol_c, unit_vec) +
                      mol_surf_rep_2DMA(ol_o, ol_c, 4)
    end

    # overlayer-monolayer interaction

   pot_mlol = Float64(0.0)
   for i in 1:nmols_ml, j in 1:nmols_ol2
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [xy_ol[j], xy_ol[j+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end
    for i in 1:nmols_ml, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], ϕ_ol[j], θ_ol[j])
    end
    return pot_mlml , pot_olol , (pot_mlsurf , pot_olsurf).*joule2wn , pot_mlol
end 

# Calculate the ith monolayer molecule contribution into the energy
function energy_ml_single(x, lattice_ml, lattice_ol, ϕ_ol, θ_ol, i) 

    theta_ml = x[1+0*nmols_ml:1*nmols_ml]
    phi_ml =   x[1+1*nmols_ml:2*nmols_ml]
    δr_ml =    x[1+2*nmols_ml:5*nmols_ml]
    δz_ol =    x[ndofs_ml]
    theta_ol = x[1+ndofs_ml+0*nmols_ol2:ndofs_ml+1*nmols_ol2]
    phi_ol   = x[1+ndofs_ml+1*nmols_ol2:ndofs_ml+2*nmols_ol2]
    xy_ol    = x[1+ndofs_ml+2*nmols_ol2:ndofs_ml+4*nmols_ol2]
    
    i_range = 1:nmols_ml

    # intra-monolayer interaction for imol

    pot_mlml = Float64(0.0)
    for j in i_range[i_range .!= i]
        rvec12 = lattice_ml[i,:] - lattice_ml[j,:] + δr_ml[i:nmols_ml:end] - δr_ml[j:nmols_ml:end] 
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlml += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end

    # Monolayer-Surface interaction for imol

    pot_mlsurf = Float64(0.0)

    stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
    rvec = [δr_ml[i], δr_ml[i+nmols_ml], δr_ml[i+2*nmols_ml] + lattice_ml[i,3]] * a0_surf
    unit_vec = [stheta*cosphi, stheta*sphi, costheta]

    ml_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    pot_mlsurf += mol_surf_attr_2DMA_tensor(ml_o, ml_c, unit_vec) +
                      mol_surf_rep_2DMA(ml_o, ml_c, 4)

    # overlayer-monolayer interaction

   pot_mlol = Float64(0.0)
   for j in 1:nmols_ol2
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [xy_ol[j], xy_ol[j+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end
    for j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], ϕ_ol[j], θ_ol[j])
    end


return pot_mlml + pot_mlol + pot_mlsurf*joule2wn
end 

# Calculate the ith overlayer molecule contribution into the energy
function energy_ol_single(x, lattice_ml, lattice_ol, ϕ_ol, θ_ol, i) 

    theta_ml = x[1+0*nmols_ml:1*nmols_ml]
    phi_ml =   x[1+1*nmols_ml:2*nmols_ml]
    δr_ml =    x[1+2*nmols_ml:5*nmols_ml]
    δz_ol =    x[ndofs_ml]
    theta_ol = x[1+ndofs_ml+0*nmols_ol2:ndofs_ml+1*nmols_ol2]
    phi_ol   = x[1+ndofs_ml+1*nmols_ol2:ndofs_ml+2*nmols_ol2]
    xy_ol    = x[1+ndofs_ml+2*nmols_ol2:ndofs_ml+4*nmols_ol2]
    
    
    i_range = 1:nmols_ol2

    # intra-overlayer interaction for imol
    # println("mol. 1","\t", "mol. 2","\t", "Distance/Å","\t \t", "Energy/cm-1")
    pot_olol = Float64(0.0)
    for j in i_range[i_range .!= i]
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] +[ xy_ol[i:nmols_ol2:end] - xy_ol[j:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_olol += co_co_interaction(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j])
    #    println(i,"\t", j,"\t", norm(rvec12/1e-10),"\t \t", co_co_interaction(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j]))
    end
    for j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] + [xy_ol[i:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        # println(rvec12)
        pot_olol += co_co_interaction(rvec12, phi_ol[i], theta_ol[i], ϕ_ol[j], θ_ol[j])
    end

    # Overlayer-Surface interaction
    pot_olsurf = Float64(0.0)
    rvec = lattice_ol[i,:]
    rvec[1] -= round(rvec[1])
    rvec[2] -= round(rvec[2])

    stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
    rvec += [xy_ol[i], xy_ol[i+nmols_ol2], δz_ol] 
    rvec *= a0_surf
    unit_vec = [stheta*cosphi, stheta*sphi, costheta]

    ol_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ol_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ol_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    pot_olsurf += mol_surf_attr_2DMA_tensor(ol_o, ol_c, unit_vec) +
                  mol_surf_rep_2DMA(ol_o, ol_c, 4)

    # overlayer-monolayer interaction
    pot_mlol = Float64(0.0)
   for j in 1:nmols_ml
        rvec12 = lattice_ml[j,:] - lattice_ol[i,:] + δr_ml[j:nmols_ml:end] - [xy_ol[i], xy_ol[i+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[j], theta_ml[j], phi_ol[i], theta_ol[i])
    end

    return pot_olol + pot_olsurf*joule2wn + pot_mlol
end 

# Calculate the energy (overlayer-monolayer + Overlayer-Surface)
function energy_ol_δz(x, lattice_ml, lattice_ol, ϕ_ol, θ_ol) 

    theta_ml = x[1+0*nmols_ml:1*nmols_ml]
    phi_ml =   x[1+1*nmols_ml:2*nmols_ml]
    δr_ml =    x[1+2*nmols_ml:5*nmols_ml]
    δz_ol =    x[ndofs_ml]
    theta_ol = x[1+ndofs_ml+0*nmols_ol2:ndofs_ml+1*nmols_ol2]
    phi_ol   = x[1+ndofs_ml+1*nmols_ol2:ndofs_ml+2*nmols_ol2]
    xy_ol    = x[1+ndofs_ml+2*nmols_ol2:ndofs_ml+4*nmols_ol2]

    # Overlayer-Surface interaction

    pot_olsurf = Float64(0.0)
    for i in 1:nmols_ol2
        
        rvec = lattice_ol[i,:]
        rvec[1] -= round(rvec[1])
        rvec[2] -= round(rvec[2])

        stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
        rvec += [xy_ol[i], xy_ol[i+nmols_ol2], δz_ol] 
        rvec *= a0_surf
        unit_vec = [stheta*cosphi, stheta*sphi, costheta]

        ol_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ol_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ol_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        pot_olsurf += mol_surf_attr_2DMA_tensor(ol_o, ol_c, unit_vec) +
                      mol_surf_rep_2DMA(ol_o, ol_c, 4)
    end

    # overlayer-monolayer interaction

   pot_mlol = Float64(0.0)
   for i in 1:nmols_ml, j in 1:nmols_ol2
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [xy_ol[j], xy_ol[j+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end
    for i in 1:nmols_ml, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], ϕ_ol[j], θ_ol[j])
    end



 return pot_mlol + pot_olsurf*joule2wn

end 
