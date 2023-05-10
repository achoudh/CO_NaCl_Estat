include("constants.jl")
#include("co_co_Guo.jl")
#include("lattice_construction.jl")
include("Plots_default.jl")
include("site_surface_interaction.jl")
using Random

function x_step(x::Vector{Float64} ,rstep::Float64, astep::Float64)

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

acceptance_probability(delta::Float64, T::Float64)::Float64 = exp(-delta/T)

function annealing_schedule(temp::Float64, cooling_rate::Float64)::Float64
    return temp * cooling_rate
end

function simulated_annealing(initial_state::Vector{Float64}, rstep::Float64, astep::Float64,  
    cooling_rate::Float64, max_temperature::Float64, max_iterations::Int64 , nift::Int64, tin::Int64)

    del::Vector{Float64} = []
    int_min = zeros(Float64,tin)
    # Initialize the current state
    current_state::Vector{Float64} = initial_state
    current_energy::Float64  = energy(current_state)
    
    # Initialize the best state
    best_state::Vector{Float64} = current_state
    best_energy::Float64  = current_energy
    
    for it::Int32 in 1:tin
    # Initialize the temperature
    temp = max_temperature
    delta::Float64 = 0.0
    
    for tk::Int64 in 1:max_iterations # Loop over the iterations

        for i::Int64 in 1:nift # number of iterations at the given temp
            # Generate a new state
            new_state::Vector{Float64} = x_step(current_state, rstep, astep)
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
    int_min[it] = best_energy 
    end
    return best_state, best_energy, del, int_min
end

# Single CO on NaCl
function energy_single(x::Vector{Float64}) ::Float64

    rvec::Vector{Float64} = x[1:3] * a0_surf
    theta::Float64 = x[4] * degrees
    phi::Float64 = x[5] * degrees

    stheta::Float64, sphi::Float64, costheta::Float64, cosphi::Float64 = sin(theta), sin(phi), cos(theta), cos(phi)

    ml_o::Vector{Float64} = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c::Vector{Float64} = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ml_bc::Vector{Float64} = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    out_attr::Float64 = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
    out_rep::Float64, out_disp::Float64  = mol_surf_rep_stone(ml_o, ml_c, 4)
    
    return (out_attr+out_rep-out_disp) * joule2wn
end
 
function energy(x) #all coordinate of the xC_unitcell is set free except z_ml

    theta_ml = x[1+0*nmols_ml:1*nmols_ml] .* pi / 180.0 #fill(38.0,4) * pi / 180.0
    phi_ml = x[1+1*nmols_ml:2*nmols_ml] .* pi / 180
    ml_in = zeros(Float64, nmols_ml, 3)
    ml_in[:,1] = x[1+2*nmols_ml:3*nmols_ml]
    ml_in[:,2] = x[1+3*nmols_ml:4*nmols_ml]
    ml_in[:,3] = x[1+4*nmols_ml:5*nmols_ml] #fill(x[17],4)
    
    ml_in[:,1] = ml_in[:,1] - round.(ml_in[:,1])
    ml_in[:,2] = ml_in[:,2] - round.(ml_in[:,2])
    
    comml = comml1 + ml_in
    
    Vmlml = zeros(Float64,10)
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        local rvec12 = (comml[i,:] - comml[j,:])
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        Vmlml .+= site_site_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i],bondlength_ml[j], phi_ml[j], theta_ml[j])
    end
    
    Vmlsurf = zeros(Float64,3)
    for i in 1:nmols_ml
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = ml_in[i,:] .* a0_surf
        rvec[3] += z_ml_ang

        ml_o = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
        out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c, 4)
        Vmlsurf .+= [out_attr, out_rep,- out_disp]
    end

    out = Vmlml[1] + sum(Vmlsurf)

   return out*joule2wn/nmols_ml
end 

initial_state::Vector{Float64} = [0.0, 0.0, 3.1e-10/a0_surf, 0.0 , 1.0] 


@time res = simulated_annealing(initial_state, 0.1, 5.0, 0.1, 10000.0, 10, 100, 3)

display(plot(res[3])) # .- res[3][1]))
println(res[1])
println(res[2])
println(res[4])

# energy(initial_state)


# using Profile , ProfileView

# @profile simulated_annealing(initial_state, 30, 0.5)


# ProfileView.view()