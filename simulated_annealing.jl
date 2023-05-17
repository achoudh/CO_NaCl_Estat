function new_coords(x::Vector{Float64}, δq::Vector{Float64}, flgs::Vector{Int32}, 
                    i::Int)

    f = flgs[i]
    x_new::Float64 = x[i] + (rand() - 0.5)*δq[i]

    if f == 1 # θ-type coordinates
        x_new = mod(x_new, 360.0)
        x_new = x_new > 180.0 ?  360.0 - x_new : x_new
    elseif f == 2 # ϕ-type coordinates
        x_new = mod(x_new, 360.0)
    elseif f == 3 # in-plane coordinates
        x_new -= round(x_new)
    end

    return x_new
end

# Method to calculate the total energy
function energy(x,lattice) 

    theta_ml = x[1+0*nmols_ml:1*nmols_ml]*degrees 
    phi_ml = x[1+1*nmols_ml:2*nmols_ml]*degrees
    ml_in =  x[1+2*nmols_ml:5*nmols_ml]
    
    pot_mlml, pot_mlsurf = Float64(0.0), Float64(0.0)
    
    # CO-CO monolayer interaction
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        rvec12 = lattice[i,:] - lattice[j,:] + ml_in[i:nmols_ml:end] - ml_in[j:nmols_ml:end] 
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlml += co_co_int_nu(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end

    # Monolayer-Surface interaction
    for i in 1:nmols_ml

        stheta, sphi, costheta, cosphi = sin(x[i]), sin(x[i+nmols_ml]), cos(x[i]), cos(x[i+nmols_ml])
        rvec = [x[i+2*nmols_ml], x[i+3*nmols_ml], x[i+4*nmols_ml] + lattice[i,3]] * a0_surf
    
        ml_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        pot_mlsurf += mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta) + 
                      mol_surf_rep_stone(ml_o, ml_c, 4)
    end

    return pot_mlml + pot_mlsurf*joule2wn
end 

# Method to calculate the ith molecule contribution into the energy
function energy(x,lattice,i) 

    theta_ml = x[1+0*nmols_ml:1*nmols_ml]*degrees 
    phi_ml = x[1+1*nmols_ml:2*nmols_ml]*degrees
    ml_in =  x[1+2*nmols_ml:5*nmols_ml]
    
    pot_mlml   = Float64(0.0)
    pot_mlsurf = Float64(0.0)
    i_range = 1:nmols_ml


    # CO-CO monolayer interaction for imol
    for j in i_range[i_range .!= i]
        rvec12 = lattice[i,:] - lattice[j,:] + ml_in[i:nmols_ml:end] - ml_in[j:nmols_ml:end] 
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlml += co_co_int_nu(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end

    # Monolayer-Surface interaction
    stheta, sphi, costheta, cosphi = sin(x[i]), sin(x[i+nmols_ml]), cos(x[i]), cos(x[i+nmols_ml])
    rvec = [x[i+2*nmols_ml], x[i+3*nmols_ml], x[i+4*nmols_ml] + lattice[i,3]] * a0_surf

    ml_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    pot_mlsurf += mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta) + 
                    mol_surf_rep_stone(ml_o, ml_c, 4)

    return pot_mlml + pot_mlsurf*joule2wn
end 

acceptance_probability(delta::Float64, T::Float64)::Float64 = exp(-delta/T)

function annealing_schedule(temp::Float64, cooling_rate::Float64)::Float64
    return temp * cooling_rate
end

function simulated_annealing(initial_state::Vector{Float64}, lattice, 
                                δq::Vector{Float64}, flgs::Vector{Int32}, 
                                cooling_rate::Float64, 
                                max_temperature::Float64, n_iterations::Int64 , 
                                nstep_thermalization::Int64, n_annealing_cycles::Int64)

    int_min = zeros(Float64,n_annealing_cycles)
    # Initialize the current state
    current_state = deepcopy(initial_state)
    current_energy::Float64  = energy(current_state,lattice)
    new_state = deepcopy(initial_state)
    del::Vector{Float64} = [current_energy]
    
    # Initialize the best state
    best_state = deepcopy(initial_state)
    best_energy::Float64  = current_energy
    delta::Float64 = 0.0

    # Exclude frozen Dofs
    flex_dofs = (1:size(initial_state,1))[flgs .!=0]

    for it in 1:n_annealing_cycles
        # Initialize the temperature
        temp = max_temperature
    
        for tk::Int64 in 1:n_iterations # Loop over the iterations

            for i::Int64 in flex_dofs # Loop over DoFs

                # Find a molecule to which ith DoF belongs
                imol = mod(i-1, nmols_ml) + 1
#println((i, imol,current_energy))
                # Get the imol-th energy contribution
                current_energy_imol::Float64  = energy(current_state, lattice, imol)
#println(current_energy_imol)
#println(current_state)
                # Generate a new state
                new_state[i] = new_coords(current_state, δq, flgs,i)
                new_energy_imol::Float64  = energy(new_state, lattice, imol)
                # Calculate the energy difference
                delta  = new_energy_imol - current_energy_imol
#println(new_energy_imol)
#println(new_state)
#println(delta)
                
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
    end
    return best_state, best_energy, del, int_min
end