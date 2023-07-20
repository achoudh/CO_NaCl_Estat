function new_coords(x::Vector{Float64}, δq::Vector{Float64}, flgs::Vector{Int32}, 
                    i::Int)

    f = flgs[i]
    x_new::Float64 = x[i] + (rand() - 0.5)*δq[i]

    if f == 1 # θ-type coordinates
        x_new = mod(x_new, 2*pi)
        x_new = x_new > pi ?  2*pi - x_new : x_new
    elseif f == 2 # ϕ-type coordinates
        x_new = mod(x_new, 2*pi)
    elseif f == 3 # in-plane coordinates
        x_new -= round(x_new)
    end

    return x_new
end

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
        pot_mlml += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end

    # intra-overlayer interaction
    println("mol. 1","\t", "mol. 2","\t", "Distance/Å","\t \t", "Energy/cm-1")
    pot_olol = Float64(0.0)
    for i in 1:nmols_ol2-1, j in i+1:nmols_ol2
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] +[ xy_ol[i:nmols_ol2:end] - xy_ol[j:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_olol += co_co_NNpes(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j])
        println(i,"\t", j,"\t",rvec12./1e-10,"\t \t", norm(rvec12)/1e-10,"\t \t", co_co_NNpes(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j]))
    end
    for i in 1:nmols_ol2, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] + [xy_ol[i:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        # println(rvec12)
        pot_olol += co_co_NNpes(rvec12, phi_ol[i], theta_ol[i], ϕ_ol[j], θ_ol[j])
    end

    # Monolayer-Surface interaction

    pot_mlsurf = Float64(0.0)
    for i in 1:nmols_ml

        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = [δr_ml[i], δr_ml[i+nmols_ml], δr_ml[i+2*nmols_ml] + lattice_ml[i,3]] * a0_surf
    
        ml_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        pot_mlsurf += mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta) + 
                      mol_surf_rep_stone(ml_o, ml_c, 4)
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

        ol_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ol_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ol_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        pot_olsurf += mol_surf_attr_stone(ol_o, ol_c, ol_bc, costheta) + 
                    mol_surf_rep_stone(ol_o, ol_c, 4)
    end

    # overlayer-monolayer interaction

   pot_mlol = Float64(0.0)
   for i in 1:nmols_ml, j in 1:nmols_ol2
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [xy_ol[j], xy_ol[j+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end
    for i in 1:nmols_ml, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], ϕ_ol[j], θ_ol[j])
    end

    return pot_mlml + pot_olol + (pot_mlsurf + pot_olsurf)*joule2wn + pot_mlol
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
        pot_mlml += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end

    # Monolayer-Surface interaction for imol

    pot_mlsurf = Float64(0.0)

    stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
    rvec = [δr_ml[i], δr_ml[i+nmols_ml], δr_ml[i+2*nmols_ml] + lattice_ml[i,3]] * a0_surf

    ml_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    pot_mlsurf += mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta) + 
                    mol_surf_rep_stone(ml_o, ml_c, 4)

    # overlayer-monolayer interaction

   pot_mlol = Float64(0.0)
   for j in 1:nmols_ol2
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [xy_ol[j], xy_ol[j+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end
    for j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], ϕ_ol[j], θ_ol[j])
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
        pot_olol += co_co_NNpes(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j])
    #    println(i,"\t", j,"\t", norm(rvec12/1e-10),"\t \t", co_co_NNpes(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j]))
    end
    for j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] + [xy_ol[i:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        println(rvec12)
        pot_olol += co_co_NNpes(rvec12, phi_ol[i], theta_ol[i], ϕ_ol[j], θ_ol[j])
    end

    # Overlayer-Surface interaction
    pot_olsurf = Float64(0.0)
    rvec = lattice_ol[i,:]
    rvec[1] -= round(rvec[1])
    rvec[2] -= round(rvec[2])

    stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
    rvec += [xy_ol[i], xy_ol[i+nmols_ol2], δz_ol] 
    rvec *= a0_surf

    ol_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ol_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ol_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    pot_olsurf += mol_surf_attr_stone(ol_o, ol_c, ol_bc, costheta) + 
                mol_surf_rep_stone(ol_o, ol_c, 4)

    # overlayer-monolayer interaction
    pot_mlol = Float64(0.0)
   for j in 1:nmols_ml
        rvec12 = lattice_ml[j,:] - lattice_ol[i,:] + δr_ml[j:nmols_ml:end] - [xy_ol[i], xy_ol[i+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_NNpes(rvec12, phi_ml[j], theta_ml[j], phi_ol[i], theta_ol[i])
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

        ol_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ol_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ol_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        pot_olsurf += mol_surf_attr_stone(ol_o, ol_c, ol_bc, costheta) + 
                    mol_surf_rep_stone(ol_o, ol_c, 4)
    end

    # overlayer-monolayer interaction

   pot_mlol = Float64(0.0)
   for i in 1:nmols_ml, j in 1:nmols_ol2
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [xy_ol[j], xy_ol[j+nmols_ol2], δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end
    for i in 1:nmols_ml, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_NNpes(rvec12, phi_ml[i], theta_ml[i], ϕ_ol[j], θ_ol[j])
    end



 return pot_mlol + pot_olsurf*joule2wn

end 

acceptance_probability(delta::Float64, T::Float64)::Float64 = exp(-delta/T)

function annealing_schedule(temp::Float64, cooling_rate::Float64)::Float64
    return temp * cooling_rate
end

function simulated_annealing(initial_state::Vector{Float64}, lattice_ml, 
                             lattice_ol, phi_ol, theta_ol, trig_uc,
                                δq::Vector{Float64}, flgs::Vector{Int32}, 
                                cooling_rate::Float64, 
                                max_temperature::Float64, n_iterations::Int64 , 
                                nstep_thermalization::Int64, n_annealing_cycles::Int64)

    int_min = zeros(Float64,n_annealing_cycles)
    # Initialize the current state
    current_state = deepcopy(initial_state)
    current_energy::Float64  = energy(current_state, lattice_ml, lattice_ol, phi_ol, theta_ol)
    new_state = deepcopy(initial_state)
    del::Vector{Float64} = [current_energy]
 
    # Initialize the best state
    best_state = deepcopy(initial_state)
    best_energy::Float64  = current_energy
    delta::Float64 = 0.0

    # Keep track of all best states
    all_best_states = zeros(Float64, n_annealing_cycles, size(initial_state,1))

    # Exclude frozen Dofs
    flex_dofs = (1:size(initial_state,1))[flgs .!=0]

    for it in 1:n_annealing_cycles
        # Initialize the temperature
        temp = max_temperature
    
        for tk::Int64 in 1:n_iterations # Loop over the iterations

            for i::Int64 in flex_dofs # Loop over DoFs

                if i < ndofs_ml # monolayer Dofs

                    # Find a molecule to which ith DoF belongs
                    imol = mod(i-1, nmols_ml) + 1
                    #println((i, imol,current_energy))
                    # Get the imol-th energy contribution
                    current_energy_imol  = 
                        energy_ml_single(current_state, lattice_ml, lattice_ol, phi_ol, theta_ol, imol)
                    #println(current_energy_imol)
                    #println(current_state)
                    # Generate a new state
                    new_state[i]         = new_coords(current_state, δq, flgs,i)
                    new_energy_imol      = 
                        energy_ml_single(new_state, lattice_ml, lattice_ol, phi_ol, theta_ol, imol)
                    # Calculate the energy difference
                    delta  = new_energy_imol - current_energy_imol
                    #println(new_energy_imol)
                    #println(new_state)
                    #println(delta)

                elseif i == ndofs_ml # Get the overlayer interaction w.r.t. δz_ol

                    # Get the δz_ol energy contribution
                    current_energy_imol  = 
                        energy_ol_δz(current_state, lattice_ml, lattice_ol, phi_ol, theta_ol)
                    # Generate a new state
                    new_state[i]         = new_coords(current_state, δq, flgs,i)
                    new_energy_imol      = 
                        energy_ol_δz(new_state, lattice_ml, lattice_ol, phi_ol, theta_ol)
                    # Calculate the energy difference
                    delta  = new_energy_imol - current_energy_imol

                else
                    # Find a molecule to which ith DoF belongs
                    imol = mod(i-1-ndofs_ml, nmols_ol2) + 1
                    #println((i, imol,current_energy))
                    # Get the imol-th energy contribution
                    current_energy_imol  = 
                        energy_ol_single(current_state, lattice_ml, lattice_ol, phi_ol, theta_ol, imol)
                    #println(current_energy_imol)
                    #println(current_state)
                    # Generate a new state
                    new_state[i]         = new_coords(current_state, δq, flgs,i)
                    new_energy_imol      = 
                        energy_ol_single(new_state, lattice_ml, lattice_ol, phi_ol, theta_ol, imol)
                    # Calculate the energy difference
                    delta  = new_energy_imol - current_energy_imol
                    #println(new_energy_imol)
                    #println(new_state)
                    #println(delta)
                end
                
                # Decide whether to accept the new state
                if delta <= 0 || rand() < acceptance_probability(delta, temp)
                    current_state[i] = new_state[i]
                    current_energy = current_energy + delta
                    push!(del, current_energy)
                else
                    new_state[i] = current_state[i]
                end
                
                # Update the best state if necessary
                if current_energy < best_energy
                    best_state[:] = current_state[:]
                    best_energy = current_energy
                end
            end

            # Update the temperature
            temp = annealing_schedule(temp, cooling_rate)
        end
        int_min[it] = best_energy 

        # Display Structure and IR Spectra
        ml_structure = structure_unitmono(best_state, com0_ml, com0_ol)
        ipda, isda, ip, is = ir_spectra(νk, best_state, com0_ml, Δν)
        ml_spectra = plot(νk, [ipda isda], label=["p-pol" "s-pol"], title=string("IR-Spectra (domain averaged)", string(it)))
        combined_plot = plot(ml_spectra, ml_structure, layout = (2, 1), size = (800, 800))
        display(combined_plot)

        all_best_states[it,:] = best_state
    end
    return best_state, best_energy, del, int_min, all_best_states
end