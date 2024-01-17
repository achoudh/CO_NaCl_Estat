using Random
using LinearAlgebra
using Printf

#########################
# Pre-defined functions #
#########################

include("constants.jl")
include("lattice_construction.jl")
#include("visualization_makie.jl")
include("write_to_file.jl")

## Turn on the interaction needed ##
include("energy_cal_new.jl") 

include("simulated_annealing.jl")
include("ir_spectra.jl")

##################
# Initialization #
##################

####################################
# Current update                   #
# Steps and initial conditions are #
# limited to only have the C-down  #
####################################


# construct a monolayer

# orientation of molecules in a monolayer's unit cell
θ_uc = zeros(Float64, 4) + [0.0,0.0,0.0,0.0]*degrees +[24.0,24.0,24.0,24.0] .* degrees
ϕ_uc = zeros(Float64, 4) + [0.0,0.0,0.0,0.0]*degrees +[360.0,360.0,180.0,180.0] .* degrees
# monolayer-surface distance (reduced units)
z_ml = 3.35e-10/a0_surf
# get a monolayer molecules' reduced positions and orientation
com0_ml, phi_ml, theta_ml = monolayer(θ_uc, ϕ_uc, z_ml)
# deviation vectors of molecular positions (r = com0 + δr)
δr_ml = zeros(Float64,nmols_ml,3)
δr_ml[1 + 0*nmols_ml:1*nmols_ml] = repeat([0.16, 0.16, -0.16, -0.16], outer =nx*ny)   # δr
δr_ml[1 + 1*nmols_ml:2*nmols_ml] = repeat([0.0, 0.0, 0.0, 0.0], outer =nx*ny)
δr_ml[1 + 2*nmols_ml:3*nmols_ml] = fill(-0.03/3.99, nmols_ml)

# construct an overlayer

# orientation of molecules in an overlayer's unit cell
θ_uc = [ 3, 1, 3, 1, 1, 3, 1, 3]*pi/4.0  # The old structure [ 3, 1, 3, 1, 3, 1, 3, 1]*pi/4.0 is not correct
ϕ_uc = [-1, 0,-1, 0, 2, 1, 2, 1]*pi/2.0 # The old structure [-1, 0,-1, 0, 1, 2, 1, 2]*pi/4.0
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
initial_state[ndofs_ml]                  = 10.0          # overlayer height deviation from c.-of-m.
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
flgs[1 + 2*nmols_ml:4*nmols_ml]   = fill(0, 2*nmols_ml)
flgs[1 + 2*nmols_ml:4:3*nmols_ml]   = fill(3, Int(nx*ny))
flgs[1 + 3*nmols_ml:4:4*nmols_ml]   = fill(3, Int(nx*ny))
# for i in 1 + 2*nmols_ml:3*nmols_ml
#     if (i-1)%4 < 2
#         flgs[i] = 30
#     else
#         flgs[i] = 31
#     end
# end


flgs[1 + 4*nmols_ml:5*nmols_ml]   = fill(4, nmols_ml) 
flgs[ndofs_ml]                    = 0 # adding the overlayer shift
flgs[1 + ndofs_ml + 0*nmols_ol2 : ndofs_ml + 1*nmols_ol2] = fill(0, nmols_ol2)
flgs[1 + ndofs_ml + 1*nmols_ol2 : ndofs_ml + 2*nmols_ol2] = fill(0, nmols_ol2)
flgs[1 + ndofs_ml + 2*nmols_ol2 : ndofs_ml + 4*nmols_ol2] = fill(0, 2*nmols_ol2)

println(energy(initial_state,com0_ml,com0_ol, phi_ol, theta_ol)  )

#######################################################
# Run simulation with randomly modified initial_state #
#######################################################
#task_id = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

#seed_num = 1234 + task_id
#Random.seed!(seed_num);

# Set step sizes

#if task_id == 1
println("initial tilt angle is 0 to pi")
println("The simulation is only run so that final structure is only C-down")
println("overlayer is far and fixed")
println("annealing parameters")
println("0.1, 10000.0, 600, 10, 10")
println("Introducing new colling mechanism and turning on the theramlization")
#end


modified_state = random_coords(initial_state,flgs,[pi/2, 2*pi, 0.5, 0.2])

# println(task_id, energy(modified_state,com0_ml,com0_ol, phi_ol, theta_ol))    

res = simulated_annealing(modified_state, com0_ml, com0_ol, phi_ol, theta_ol, trig_uc, 
                            δq, flgs, 
                            0.4, 1000.0, 1, 1, 1) ## The thermalization loop does not matter


# write_to_file("./Output/22-08-2023/2/buried_ov_fixed_dof_$task_id.txt", res)



res[2][end]