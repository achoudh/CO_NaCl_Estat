################################################################################
                                    ########
                                    # NOTE #
                                    ########

# This program finds the minimum energy structure for bare and buried monolayer
# Only considers the attraction interactions 
# includes "attraction_PES_arnab.jl" where the attraction protections are listed


################################################################################

using Random
using LinearAlgebra
using Printf
using Optim

#########################
# Pre-defined functions #
#########################

include("constants.jl")
include("lattice_construction.jl")
include("visualization_makie.jl")
include("ir_spectra.jl")
include("attraction_PES_arnab.jl")


##################
# Initialization #
##################

# Random.seed!(1234);

# construct a monolayer

# orientation of molecules in a monolayer's unit cell
θ_uc = zeros(Float64, 4) + [30.0,30.0,30.0,30.0] * degrees
ϕ_uc = zeros(Float64, 4) + [20.0,60.0,20.0,60.0] * degrees
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

###########################
# set up initial geometry #
###########################

# Set initial geometry for monolayer
initial_state = zeros(Float64, 2*nmols_ml)
initial_state[1 + 0*nmols_ml:1*nmols_ml] = theta_ml     # θ
initial_state[1 + 1*nmols_ml:2*nmols_ml] = phi_ml       # ϕ 

δr_ml = vec(δr_ml)   # δr
# δr_ml[1 + 0*nmols_ml:1*nmols_ml] = vec(δr_ml)   # δr
# δr_ml[1 + 1*nmols_ml:2*nmols_ml] = vec(δr_ml)   # δr
# δr_ml[1 + 2*nmols_ml:3*nmols_ml] = fill(δz_ml, nmols_ml)   # δr

lower = zeros(Float64, 2*nmols_ml)
upper = zeros(Float64, 2*nmols_ml)
upper[1 + 0*nmols_ml:1*nmols_ml] = fill(π, nmols_ml)     # θ
upper[1 + 1*nmols_ml:2*nmols_ml] = fill(2*π, nmols_ml)


# overlayer
overlayer_states = zeros(Float64, 2*nmols_ol2)
overlayer_states[1 + 0*nmols_ol2 : 1*nmols_ol2] = theta_ol[1:nmols_ol2]    # θ
overlayer_states[1 + 1*nmols_ol2 : 2*nmols_ol2] = phi_ol[1:nmols_ol2]    # ϕ 

###################
# Energy function #
###################

# Total energy (intra-monolayer + intra-overlayer + monolayer-Surface + overlayer-Surface + monolayer-overlayer )
function energy(x::Vector{Float64}, lattice_ml::Matrix{Float64}=com0_ml,
             lattice_ol::Matrix{Float64}=com0_ol, ϕ_ol::Vector{Float64}=phi_ol, θ_ol::Vector{Float64}=theta_ol) 

    theta_ml = x[1+0*nmols_ml:1*nmols_ml]
    phi_ml   = x[1+1*nmols_ml:2*nmols_ml]

    δr_ml[1 + 2*nmols_ml:3*nmols_ml] = fill(δz_ml, nmols_ml)
    
    theta_ol = overlayer_states[1+0*nmols_ol2:1*nmols_ol2]
    phi_ol   = overlayer_states[1+1*nmols_ol2:2*nmols_ol2]
    # xy_ol    = fixed_states[1+2*nmols_ol2:4*nmols_ol2]
    
    # intra-monolayer interaction

    pot_mlml = Float64(0.0)
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        rvec12 = lattice_ml[i,:] - lattice_ml[j,:] + δr_ml[i:nmols_ml:end] - δr_ml[j:nmols_ml:end] 
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlml += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end

    # intra-overlayer interaction
    # println("mol. 1","\t", "mol. 2","\t", "Distance/Å","\t \t", "Energy/cm-1")
    pot_olol = Float64(0.0)
    for i in 1:nmols_ol2-1, j in i+1:nmols_ol2
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] # +[ xy_ol[i:nmols_ol2:end] - xy_ol[j:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_olol += co_co_interaction(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j])
        # println(i,"\t", j,"\t", norm(rvec12)/1e-10,"\t \t", co_co_interaction(rvec12, phi_ol[i], theta_ol[i], phi_ol[j], theta_ol[j]))
    end
    # println(pot_olol)
    for i in 1:nmols_ol2, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ol[i,:] - lattice_ol[j,:] # + [xy_ol[i:nmols_ol2:end] ; 0 ]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        # println(rvec12)
        pot_olol += co_co_interaction(rvec12, phi_ol[i], theta_ol[i], ϕ_ol[j], θ_ol[j])
    end

    # Monolayer-Surface interaction

    pot_mlsurf = Float64(0.0)
    for i in 1:nmols_ml

        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = [δr_ml[i], δr_ml[i+nmols_ml], δr_ml[i+2*nmols_ml] + lattice_ml[i,3]] * a0_surf
                
        pot_mlsurf += mol_surf_attr_arnab(rvec, costheta)
    end

    # Overlayer-Surface interaction

    pot_olsurf = Float64(0.0)
    for i in 1:nmols_ol2

        rvec = lattice_ol[i,:]
        rvec[1] -= round(rvec[1])
        rvec[2] -= round(rvec[2])

        stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
        rvec += [0.0, 0.0, δz_ol] 
        rvec *= a0_surf

        pot_olsurf += mol_surf_attr_arnab(rvec, costheta)
    end

    # overlayer-monolayer interaction

   pot_mlol = Float64(0.0)
   for i in 1:nmols_ml, j in 1:nmols_ol2
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end
    for i in 1:nmols_ml, j in 1+nmols_ol2:nmols_ol
        rvec12 = lattice_ml[i,:] - lattice_ol[j,:] + δr_ml[i:nmols_ml:end] - [0.0, 0.0, δz_ol]
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        pot_mlol += co_co_interaction(rvec12, phi_ml[i], theta_ml[i], ϕ_ol[j], θ_ol[j])
    end

    return pot_mlml + pot_olol + (pot_mlsurf + pot_olsurf)*joule2wn + pot_mlol
end 


####################
# Run optimization #
####################

δz_ml, δz_ol = 0.0, 0.0
println(energy(initial_state))

g_tol = 1e-8
x_tol = 1e-8
f_tol = 1e-8
inner_optimizer = LBFGS()
@time res_over = optimize(energy, lower, upper, initial_state, Fminbox(inner_optimizer), Optim.Options(g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, iterations = 2000))
println(res_over)
# println(res_over.minimizer)
# energy(res_over.minimizer)


# Set combined geometry
final_state = zeros(Float64, ndofs_ml+4*nmols_ol2)
# monolayer
final_state[1 + 0*nmols_ml:2*nmols_ml] = res_over.minimizer     # θ
# final_state[1 + 1*nmols_ml:2*nmols_ml] = phi_ml       # ϕ 
final_state[1 + 2*nmols_ml:5*nmols_ml] = vec(δr_ml)   # δr
final_state[ndofs_ml]                  = 0.0          # overlayer height deviation from c.-of-m.
# overlayer
final_state[1 + ndofs_ml + 0*nmols_ol2 : ndofs_ml + 1*nmols_ol2] = theta_ol[1:nmols_ol2]    # θ
final_state[1 + ndofs_ml + 1*nmols_ol2 : ndofs_ml + 2*nmols_ol2] = phi_ol[1:nmols_ol2]    # ϕ 
final_state[1 + ndofs_ml + 2*nmols_ol2 : ndofs_ml + 4*nmols_ol2] = vec(δr_ol)   # δr


show_params(final_state)# Display final Structure and IR Spectra
ipda, isda, ip, is = ir_spectra(νk, final_state, com0_ml, Δν)

fig = Figure()
ax = Axis
    (fig[1, 1],
        title = "A Makie Axis",
        xlabel = "The x label",
        ylabel = "The y label"
    )

ax = lines(νk, ipda)
ml_spectra = lines!(νk, isda)

ml_structure       = structure_unitmono(final_state, com0_ml, com0_ol)
ml_structure1      = scatter3d(ml_structure, camera=(10,20,), label = nothing, ticks=nothing, axes=nothing,zlims =(0,7))
combined_plot1     = scatter3d!(ml_structure,ml_structure1,layout=(1,2))
combined_plot      = plot!(ml_spectra, combined_plot1, layout=(2,1), size = (700, 400), dpi = 1200)

display(combined_plot)

νk, ipda
f