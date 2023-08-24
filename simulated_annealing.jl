#################################
# Method of energy calculations #
# are moved to another file     #
#################################


function new_coords(x::Vector{Float64}, δq::Vector{Float64}, flgs::Vector{Int32}, 
                    i::Int)

    f = flgs[i]
    x_new::Float64 = x[i] + (rand() - 0.5)*δq[i]

    if f == 1 # θ-type coordinates
        x_new = mod(x_new, 2*pi)
        x_new = x_new > pi ?  2*pi - x_new : x_new
        x_new = x_new > pi/2 ?  pi - x_new : x_new  
    elseif f == 2 # ϕ-type coordinates
        x_new = mod(x_new, 2*pi)
    elseif f == 3 # in-plane coordinates
        x_new -= round(x_new)
    end

    return x_new
end

function random_coords(x::Vector{Float64}, flgs::Vector{Int32}, δq)

    modified_state = zeros(Float64, ndofs_ml+4*nmols_ol2)

    for (i, f) in enumerate(flgs)
        r = rand()
        if f == 0 # Frozen
            x_new::Float64 = x[i]
        elseif f == 1 # θ-type coordinates
            x_new = x[i] + δq[1]*r # As it should be 0 to 180
        elseif f == 2 # ϕ-type coordinates
            x_new = x[i] + δq[2]*r # As it should be 0 to 360
        elseif f == 3
            x_new = x[i] + (r - δq[3]) # As it should be -0.5 to 0.5
        elseif f == 4  # z coordinate
            x_new = x[i] + (r - δq[4])  # How much should it move
        else 
            error("random_coords: Unknown flag")
        end
        modified_state[i] = x_new
    end

    return modified_state
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

    
    # Initialize the current state
    current_state = deepcopy(initial_state)
    current_energy::Float64  = energy(current_state, lattice_ml, lattice_ol, phi_ol, theta_ol)
    new_state = deepcopy(initial_state)

    int_min::Vector{Float64} = [current_energy] # Intermediate minimums or best energies
    accepted_step_energies::Vector{Float64} = [current_energy] # Energies at the accepted states
 
    # Initialize the best state
    best_state = deepcopy(initial_state)
    best_energy::Float64  = current_energy
    delta::Float64 = 0.0

    # Keep track of all best states
    all_best_states::Vector{Vector{Float64}} = [initial_state]

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
                    push!(accepted_step_energies, current_energy)
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

        push!(int_min, best_energy)
        push!(all_best_states, best_state)
    end
    return int_min, all_best_states, accepted_step_energies # Removed "best_state, best_energy" as the int_min and all_best_states already includes initial_state and all
end