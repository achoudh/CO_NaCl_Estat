using Random
using LinearAlgebra
using DelimitedFiles

include("constants.jl")
include("co_co_Guo.jl")
include("lattice_construction.jl")
include("Plots_default.jl")
include("site_surface_interaction.jl")
include("simulated_annealing.jl")



# initialize monolayer DoFs and flags showing frozen DoFs
initial_state = zeros(Float64, 5*nmols_ml)
frozen_flgs = fill(true, 5*nmols_ml)

theta_m = zeros(Float64, 4)
phi_m = zeros(Float64, 4)
z_ml = 3.1e-10/a0_surf
com0_ml, bondlength_ml, phi_ml, theta_ml = monolayer(theta_m, phi_m, z_ml)
δr_ml = zeros(Float64,nmols_ml,3)

initial_state[1 + 0*nmols_ml:1*nmols_ml] = theta_ml
initial_state[1 + 1*nmols_ml:2*nmols_ml] = phi_ml
initial_state[1 + 2*nmols_ml:5*nmols_ml] = vec(δr_ml)

@time res = simulated_annealing(initial_state, frozen_flgs, com0_ml, 
                                0.1, 5.0, 0.1, 10000.0, 10, 100, 3)

display(plot(res[3])) # .- res[3][1]))
println(res[1])
println(res[2])
println(res[4])

# energy(initial_state)


# using Profile , ProfileView

# @profile simulated_annealing(initial_state, 30, 0.5)


# ProfileView.view()
