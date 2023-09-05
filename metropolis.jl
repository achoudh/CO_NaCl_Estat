using Random
using LinearAlgebra
using Printf
using GLMakie

# GLMakie.activate!()
#########################
# Pre-defined functions #
#########################

include("constants.jl")
include("lattice_construction.jl")
include("visualization_makie.jl")

## Turn on the interaction needed ##
# include("energy_cal_hoang.jl") 
include("energy_cal_new.jl")

include("simulated_annealing.jl")
include("ir_spectra.jl")

##################
# Initialization #
##################

Random.seed!(1234);

# construct a monolayer

# orientation of molecules in a monolayer's unit cell
θ_uc = zeros(Float64, 4) + [0.0,0.0,0.0,0.0]*degrees #[30.0,30.0,30.0,30.0]
ϕ_uc = zeros(Float64, 4) + [0.0,0.0,0.0,0.0]*degrees #[20.0,60.0,20.0,60.0]
# monolayer-surface distance (reduced units)
z_ml = 3.35e-10/a0_surf
# get a monolayer molecules' reduced positions and orientation
com0_ml, phi_ml, theta_ml = monolayer(θ_uc, ϕ_uc, z_ml)
# deviation vectors of molecular positions (r = com0 + δr)
δr_ml = zeros(Float64,nmols_ml,3)

# construct an overlayer

# orientation of molecules in an overlayer's unit cell
θ_uc = [ 3, 1, 3, 1, 1, 3, 1, 3]*pi/4.0  # The old structure [ 3, 1, 3, 1, 3, 1, 3, 1]*pi/4.0 is not correct
ϕ_uc = [-1, 0,-1, 0, 1, 2, 1, 2]*pi/2.0 
trig_uc = (sin.(θ_uc), cos.(θ_uc), sin.(ϕ_uc), cos.(ϕ_uc))
# overlayer-surface distance (reduced units)
z_ol = z_ml + 0.5*a0_CO/a0_surf #+ 10.00*a0_CO/a0_surf
# get an overlayer molecules' reduced positions and orientation
com0_ol, phi_ol, theta_ol = overlayer(θ_uc, ϕ_uc, z_ol)
# deviation vectors of molecular positions (r = com0 + δr)
δr_ol = zeros(Float64,nmols_ol2,2)

###########################################
# step vector, flags and initial geometry #
###########################################

# step vector, and flags defining frozen DoFs
δq     = zeros(Float64, 5*nmols_ml+1+4*nmols_ol2)
flgs   = zeros(Int32, 5*nmols_ml+1+4*nmols_ol2)

# Set initial geometry
initial_state = zeros(Float64, ndofs_ml+4*nmols_ol2)
# monolayer
initial_state[1 + 0*nmols_ml:1*nmols_ml] = theta_ml     # θ
initial_state[1 + 1*nmols_ml:2*nmols_ml] = phi_ml       # ϕ 
initial_state[1 + 2*nmols_ml:5*nmols_ml] = vec(δr_ml)   # δr
initial_state[ndofs_ml]                  = 0.0          # overlayer height deviation from c.-of-m.
# overlayer
initial_state[1 + ndofs_ml + 0*nmols_ol2 : ndofs_ml + 1*nmols_ol2] = theta_ol[1:nmols_ol2]    # θ
initial_state[1 + ndofs_ml + 1*nmols_ol2 : ndofs_ml + 2*nmols_ol2] = phi_ol[1:nmols_ol2]    # ϕ 
initial_state[1 + ndofs_ml + 2*nmols_ol2 : ndofs_ml + 4*nmols_ol2] = vec(δr_ol)   # δr

# Set step sizes
δθ = 10.0*degrees
δϕ = 10.0*degrees
δxy = 0.05      # surface lattice units
δz = 0.1        # surface lattice units
δz_ol = 0.1     # surface lattice units
δq[1 + 0*nmols_ml:1*nmols_ml] = fill(δθ, nmols_ml)
δq[1 + 1*nmols_ml:2*nmols_ml] = fill(δϕ, nmols_ml)
δq[1 + 2*nmols_ml:4*nmols_ml] = fill(δxy, 2*nmols_ml)
δq[1 + 4*nmols_ml:5*nmols_ml] = fill(δz, nmols_ml)
δq[ndofs_ml] = δz_ol
δq[1 + ndofs_ml + 0*nmols_ol2 : ndofs_ml + 1*nmols_ol2] = fill(δθ, nmols_ol2)
δq[1 + ndofs_ml + 1*nmols_ol2 : ndofs_ml + 2*nmols_ol2] = fill(δϕ, nmols_ol2)
δq[1 + ndofs_ml + 2*nmols_ol2 : ndofs_ml + 4*nmols_ol2] = fill(δxy, 2*nmols_ol2)

# Set coordinate type flags:
# 0 means that  a frozen coordinate
# 1             a θ-coordinate
# 2             a ϕ-coordinate
# 3             a in-plane coordinate
# 4             a z-coordinate
flgs[1 + 0*nmols_ml:1*nmols_ml]   = fill(1, nmols_ml)
flgs[1 + 1*nmols_ml:2*nmols_ml]   = fill(2, nmols_ml)
flgs[1 + 2*nmols_ml:4*nmols_ml]   = fill(3, 2*nmols_ml)
flgs[1 + 4*nmols_ml:5*nmols_ml]   = fill(4, nmols_ml) 
flgs[ndofs_ml]                    = 4 # adding the overlayer shift
flgs[1 + ndofs_ml + 0*nmols_ol2 : ndofs_ml + 1*nmols_ol2] = fill(0, nmols_ol2)
flgs[1 + ndofs_ml + 1*nmols_ol2 : ndofs_ml + 2*nmols_ol2] = fill(0, nmols_ol2)
flgs[1 + ndofs_ml + 2*nmols_ol2 : ndofs_ml + 4*nmols_ol2] = fill(0, 2*nmols_ol2)


println("Initial state:")
#println(initial_state)
println(energy(initial_state,com0_ml,com0_ol, phi_ol, theta_ol))

arnab
# Display Structure and IR Spectra
fig = show_figure(initial_state, com0_ml, com0_ol, "Initial")
# save("C:/Users/achoudh/ownCloud/my work/CO_NaCl-estat/Estat_results/initial.png", fig)
# display(fig)


#######################################################
# Run simulation with randomly modified initial_state #
#######################################################


# Set step sizes

@time res = simulated_annealing(initial_state, com0_ml, com0_ol, phi_ol, theta_ol, trig_uc, 
                                δq, flgs, 
                                0.4, 1000000.0, 100, 1, 2)
                                # (initial_state::Vector{Float64}, lattice_ml, lattice_ol, phi_ol, theta_ol, trig_uc,
                                # δq::Vector{Float64}, flgs::Vector{Int32}, 
                                # cooling_rate::Float64, 
                                # max_temperature::Float64, n_iterations::Int64 , 
                                # nstep_thermalization::Int64, n_annealing_cycles::Int64)

# fig2 = show_figure(res[1], com0_ml, com0_ol, "Final")
# save("C:/Users/achoudh/ownCloud/my work/CO_NaCl-estat/Estat_results/final.png", fig2)
# display(GLMakie.Screen(), fig2)

# modified_states = []

arnab
# write_to_file("buried_ov_fixed_dof.txt", res)

Threads.@threads for i in 1:4

    println(Threads.threadid())

    modified_state = random_coords(initial_state,flgs,[π, 2*π, 0.5, 0.2])
    # push!(modified_state, modified_states)
    # fig = show_figure(modified_state, com0_ml, com0_ol, "Ininal$i")
    # save("C:/Users/achoudh/ownCloud/my work/CO_NaCl-estat/Estat_results/initial$i.png", show_figure(modified_state, com0_ml, com0_ol, "Ininal$i"))

    res = simulated_annealing(modified_state, com0_ml, com0_ol, phi_ol, theta_ol, trig_uc, 
                                δq, flgs, 
                                0.4, 1000000.0, 100, 1, 2)

    println("file$i.txt")

    # fig = show_figure(res[1], com0_ml, com0_ol, "Final$i")
    # display(GLMakie.Screen(), fig)
    # save("C:/Users/achoudh/ownCloud/my work/CO_NaCl-estat/Estat_results/Final$i.png", fig)
end


# using Profile , ProfileView

# @profile simulated_annealing(initial_state, 30, 0.5)


# ProfileView.view()
