include("constants.jl")
#include("co_co_Guo.jl")
#include("lattice_construction.jl")
include("Plots_default.jl")
include("site_surface_interaction.jl")
using Random

function x_step(x::Vector{Float64})
    rstep::Float64, astep::Float64 = 0.2, 10.0
    x_f = zeros(Float64, length(x))
    
    for i in 1:2
        x_f[i] = x[i] + (rand() - 0.5) * rstep
        x_f[i] = x_f[i] - round(x_f[i])
    end

    x_f[3] = x[3] + (rand() - 0.5) * rstep
    for i in 4:5
        x_f[i] = x[i] + (rand() - 0.5) * astep
        x_f[i] = mod(x_f[i], 360.0)
    end

    x_f[4] = x_f[4] > 180 ?  360 - x_f[4] : x_f[4]

    return x_f    
end

function ang_range(x, min_val, max_val, astep)
    x_f = x + (rand() - 0.5) * astep
    while x_f < min_val || x_f > max_val  # check if x is within range
        x_f = x + (rand() - 0.5) * astep # generate a new random value for x if it's not within range
    end   
    return x_f
end

acceptance_probability(delta::Float64, T::Float64)::Float64 = exp(-delta/T)

function annealing_schedule(temp::Float64, cooling_rate::Float64)::Float64
    return temp * cooling_rate
end

function simulated_annealing(initial_state::Vector{Float64}, max_iterations::Int64, cooling_rate::Float64)
    del::Vector{Float64} = []
    int_min = zeros(Float64,3)
    # Initialize the current state
    current_state::Vector{Float64} = initial_state
    current_energy::Float64  = energy(current_state)
    
    # Initialize the best state
    best_state::Vector{Float64} = current_state
    best_energy::Float64  = current_energy
    
    for tin in 1:3
    # Initialize the temperature
    temp::Float64 = 100000.0
    delta::Float64 = 0.0
    
    for tk in 1:max_iterations

    # Loop over the iterations
        for i::Float64 in 1:100 # max_iterations
            # Generate a new state
            new_state::Vector{Float64} = x_step(current_state)
            new_energy::Float64  = energy(new_state)
            
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
    int_min[tin] = best_energy 
    end
    return best_state, best_energy, del, int_min
end


function energy(x::Vector{Float64}) ::Float64
    rvec::Vector{Float64} = x[1:3] * a0_surf
    theta::Float64 = x[4] * pi / 180
    phi::Float64 = x[5] * pi / 180

    rvec[3] += 3.1e-10 
    stheta::Float64, sphi::Float64, costheta::Float64, cosphi::Float64 = sin(theta), sin(phi), cos(theta), cos(phi)

    ml_o::Vector{Float64} = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c::Vector{Float64} = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ml_bc::Vector{Float64} = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    out_attr::Float64 = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
    out_rep::Float64, out_disp::Float64  = mol_surf_rep_stone(ml_o, ml_c, 4)
    
    return (out_attr+out_rep-out_disp) * joule2wn
end

initial_state::Vector{Float64} = [0.0,0.0,0.0,180.0,1.0] 


@time res = simulated_annealing(initial_state, 300, 0.2)

display(plot(res[3])) # .- res[3][1]))
println(res[1])
println(res[2])
println(res[4])

# energy(initial_state)


# using Profile , ProfileView

# @profile simulated_annealing(initial_state, 30, 0.5)


# ProfileView.view()