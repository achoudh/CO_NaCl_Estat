using Random
using LinearAlgebra
using DelimitedFiles

include("constants.jl")
include("co_co_Guo.jl")
include("lattice_construction.jl")
include("Plots_default.jl")
include("site_surface_interaction.jl")
include("simulated_annealing.jl")

# Initialization
Random.seed!(1234);

# monolayer DoFs, step vector, and flags defining frozen DoFs
initial_state = zeros(Float64, 5*nmols_ml)
δq = zeros(Float64, 5*nmols_ml)
flgs = zeros(Int32, 5*nmols_ml)

# construct a monolayer

# orientation of molecules in a monolayer's unit cell
θ_uc = zeros(Float64, 4) + fill(30,4)*degrees
ϕ_uc = zeros(Float64, 4) + [0,180,0,180]*degrees
# monolayer-surface distance (reduced units)
z_ml = 3.35e-10/a0_surf
# get a monolayer molecules' reduced positions and orientation
com0_ml, phi_ml, theta_ml = monolayer(θ_uc, ϕ_uc, z_ml)
# deviation vectors of molecular positions (r = com0 + δr)
δr_ml = zeros(Float64,nmols_ml,3)

# construct an overlayer

# orientation of molecules in an overlayer's unit cell
θ_uc = [ 3, 1, 3, 1, 3, 1, 3, 1]*pi/4.0
ϕ_uc = [-1, 0,-1, 0, 1, 2, 1, 2]*pi/2.0
# monolayer-overlayer distance (in units of CO lattice constant)
dz = 0.5
# get an overlayer molecules' reduced positions and orientation
com_ol, phi_ol, theta_ol = overlayer(θ_uc, ϕ_uc, z_ml, dz)

# Set initial geometry

initial_state[1 + 0*nmols_ml:1*nmols_ml] = theta_ml     # θ
initial_state[1 + 1*nmols_ml:2*nmols_ml] = phi_ml       # ϕ 
initial_state[1 + 2*nmols_ml:5*nmols_ml] = vec(δr_ml)   # δr

# Set step sizes
δθ = 10.0*degrees
δϕ = 10.0*degrees
δxy = 0.05   # surface lattice units
δz = 0.1    # surface lattice units
δq[1 + 0*nmols_ml:1*nmols_ml] = fill(δθ, nmols_ml)
δq[1 + 1*nmols_ml:2*nmols_ml] = fill(δϕ, nmols_ml)
δq[1 + 2*nmols_ml:4*nmols_ml] = fill(δxy, 2*nmols_ml)
δq[1 + 4*nmols_ml:5*nmols_ml] = fill(δz, nmols_ml)

# Set coordinate type flags:
# 0 means that  a frozen coordinate
# 1             a θ-coordinate
# 2             a ϕ-coordinate
# 3             a in-plane coordinate
# 4             a z-coordinate
flgs[1 + 0*nmols_ml:1*nmols_ml] = fill(1, nmols_ml)
flgs[1 + 1*nmols_ml:2*nmols_ml] = fill(2, nmols_ml)
flgs[1 + 2*nmols_ml:4*nmols_ml] = fill(3, 2*nmols_ml)
flgs[1 + 4*nmols_ml:5*nmols_ml] = fill(4, nmols_ml)

println("Initial state:")
#println(initial_state)
println(energy(initial_state,com0_ml))

@time res = simulated_annealing(initial_state, com0_ml, δq, flgs, 
                            0.1, 100000.0, 200, 1, 3)

display(plot(res[3])) # .- res[3][1]))
println(res[1])
println(res[2])
println(res[4])

# energy(initial_state)


# using Profile , ProfileView

# @profile simulated_annealing(initial_state, 30, 0.5)


# ProfileView.view()
