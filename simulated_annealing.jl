function new_coords(x::Vector{Float64}, δq::Vector{Float64}, flgs::Vector{Int32})

    x_new = zeros(Float64, size(x,1))
    
    for (i,f) in enumerate(flgs)
        if f==0 # a frozen coordinate
            x_new[i] = x[i]
        else
            x_new[i] = x[i] + (rand() - 0.5)*δq[i]
            if f == 1 # θ-type coordinates
                x_new[i] = mod(x_new[i], 360.0)
                x_new[i] = x_new[i] > 180.0 ?  360.0 - x_new[i] : x_new[i]
            elseif f == 2 # ϕ-type coordinates
                x_new[i] = mod(x_new[i], 360.0)
            elseif f == 3 # in-plane coordinates
                x_new[i] -= round(x_new[i])
            end
        end
    end

    return x_new
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

    del::Vector{Float64} = []
    int_min = zeros(Float64,n_annealing_cycles)
    # Initialize the current state
    current_state::Vector{Float64} = initial_state
    current_energy::Float64  = energy(current_state,lattice)
    
    # Initialize the best state
    best_state::Vector{Float64} = current_state
    best_energy::Float64  = current_energy
    
    for it::Int32 in 1:n_annealing_cycles
    # Initialize the temperature
    temp = max_temperature
    delta::Float64 = 0.0
    
        for tk::Int64 in 1:n_iterations # Loop over the iterations

            for i::Int64 in 1:nstep_thermalization # number of iterations at the given temp
                # Generate a new state
                new_state::Vector{Float64} = new_coords(current_state, δq, flgs)
                new_energy::Float64  = energy(new_state, lattice)
#println((current_energy, new_energy))

                # Calculate the energy difference
                delta  = new_energy - current_energy

                # Decide whether to accept the new state
                if delta <= 0 || rand() < acceptance_probability(delta, temp)
                    current_state = new_state
                    current_energy = new_energy
                    push!(del, current_energy)
                end
                
                # Update the best state if necessary
                if current_energy < best_energy
                    best_state = current_state
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
