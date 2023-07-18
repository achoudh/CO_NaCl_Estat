using Random
using LinearAlgebra
using DelimitedFiles

include("constants.jl")
include("co_co_Guo.jl")
include("lattice_construction.jl")
include("visualization.jl")
include("site_surface_interaction.jl")
include("simulated_annealing.jl")
include("ir_spectra.jl")

# Initialization
Random.seed!(1234);

# monolayer DoFs, step vector, and flags defining frozen DoFs
initial_state = zeros(Float64, 5*nmols_ml+1)
δq = zeros(Float64, 5*nmols_ml+1)
flgs = zeros(Int32, 5*nmols_ml+1)

# construct a monolayer

# orientation of molecules in a monolayer's unit cell
θ_uc = zeros(Float64, 4) + [30.0,30.0,30.0,30.0]*degrees
ϕ_uc = zeros(Float64, 4) + [20.0,60.0,20.0,60.0]*degrees
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
trig_uc = (sin.(θ_uc), cos.(θ_uc), sin.(ϕ_uc), cos.(ϕ_uc))
# overlayer-surface distance (reduced units)
z_ol = z_ml + 0.5*a0_CO/a0_surf #+ 10.00*a0_CO/a0_surf
# get an overlayer molecules' reduced positions and orientation
com0_ol, phi_ol, theta_ol = overlayer(θ_uc, ϕ_uc, z_ol)

# Set initial geometry

initial_state[1 + 0*nmols_ml:1*nmols_ml] = theta_ml     # θ
initial_state[1 + 1*nmols_ml:2*nmols_ml] = phi_ml       # ϕ 
initial_state[1 + 2*nmols_ml:5*nmols_ml] = vec(δr_ml)   # δr
initial_state[1 + 5*nmols_ml]            = 0.0          # overlayer height deviation from c.-of-m.

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
δq[1 + 5*nmols_ml] = δz_ol

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
flgs[1 + 5*nmols_ml]              = 0 # adding the overlayer shift

println("Initial state:")
#println(initial_state)
println(energy(initial_state,com0_ml,com0_ol, phi_ol, theta_ol, trig_uc))

# Display Structure and IR Spectra
ml_structure = structure_unitmono(initial_state, com0_ml)
ipda, isda, ip, is = ir_spectra(νk, initial_state, com0_ml, Δν)
ml_spectra = plot(νk, [ipda isda], label=["p-pol" "s-pol"],xlabel = "Frequnecy/cm-1",title="IR-Spectra (domain averaged) initial state")
combined_plot = plot(ml_spectra, ml_structure, layout = (2, 1), size = (800, 800))
display(combined_plot)


@time res = simulated_annealing(initial_state, com0_ml,com0_ol, phi_ol, theta_ol, trig_uc, 
                                δq, flgs, 
                            0.9, 100000.0, 200, 200, 3)
                            # (initial_state::Vector{Float64}, lattice_ml, lattice_ol, phi_ol, theta_ol, trig_uc,
                            # δq::Vector{Float64}, flgs::Vector{Int32}, 
                            # cooling_rate::Float64, 
                            # max_temperature::Float64, n_iterations::Int64 , 
                            # nstep_thermalization::Int64, n_annealing_cycles::Int64)

display(plot(res[3], color = :black, label = " ", xlabel = "accepted steps", ylabel = "energy/cm-1")) # .- res[3][1]))
println(res[1])
println(res[2])
println(res[4])

# Display final Structure and IR Spectra
ml_structure = structure_unitmono(res[1], com0_ml)
ipda, isda, ip, is = ir_spectra(νk, res[1], com0_ml, Δν)
ml_spectra = plot(νk, [ipda isda], label=["p-pol" "s-pol"],xlabel = "Frequnecy/cm-1",title="IR-Spectra (domain averaged) final state")
combined_plot = plot(ml_spectra, ml_structure, layout = (2, 1), size = (800, 800))
display(combined_plot)
# energy(initial_state)

println("The Final ML parameters are")
show_params(res[1])

# using Profile , ProfileView

# @profile simulated_annealing(initial_state, 30, 0.5)


# ProfileView.view()
