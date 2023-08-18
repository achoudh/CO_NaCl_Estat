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
using GLMakie

#########################
# Pre-defined functions #
#########################

include("../constants.jl")
include("../lattice_construction.jl")
include("../visualization_makie.jl")
include("../ir_spectra.jl")
include("attraction_PES_arnab.jl")

filepath = "C:/Users/achoudh/ownCloud/my work/CO_NaCl-estat/Estat_results"

#############
# Functions #
#############


function low_high_limits(x::Vector{Float64}, flgs::Vector{Int64})
    low  = zeros(Float64, size(x, 1) )
    high = zeros(Float64, size(x, 1) )

    for (i, f) in enumerate(flgs)
        if f == 1 # θ-type coordinates
            low[i] = 0
            high[i] = π
        elseif f == 2 # ϕ-type coordinates
            low[i] = 0
            high[i] = 2*π
        elseif f == 3
            low[i] = -0.5
            high[i] = 0.5
        elseif f == 4  # z coordinate
            low[i] = -0.2
            high[i] = 0.8
        else 
            error("random_coords: Unknown flag")
        end
    end
    return low, high
end

function random_coords(x::Vector{Float64}, flgs::Vector{Int64}, δq)

    modified_state = zeros(Float64, size(x,1) )

    for (i, f) in enumerate(flgs)
        r = rand()
        if f == 0 # Frozen
            x_new::Float64 = x[i]
        elseif f == 1 # θ-type coordinates
            x_new = δq[1]*r # As it should be 0 to 180
        elseif f == 2 # ϕ-type coordinates
            x_new = δq[2]*r # As it should be 0 to 360
        elseif f == 3
            x_new = (r - δq[3]) # As it should be -0.5 to 0.5
        elseif f == 4  # z coordinate
            x_new = (r - δq[4])  # How much should it move
        else 
            error("random_coords: Unknown flag")
        end
        modified_state[i] = x_new
    end

    return modified_state
end

function optimize_function(initial_state::Vector{Float64}, lower::Vector{Float64}, upper::Vector{Float64}, energy_function::Function)

    g_tol = 1e-8
    x_tol = 1e-8
    f_tol = 1e-8
    
    inner_optimizer = LBFGS()
    result = optimize(energy_function, lower, upper, initial_state, Fminbox(inner_optimizer), 
            Optim.Options(g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, iterations = 2000))

    return result
    
end

# Calculates the vector distance between two state vectors 
vec_norm(x_new, x) = norm(x_new - x) / size(x, 1)

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


##################
# Initialization #
##################

    Random.seed!(1234);

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
    δr_ml[1 + 0*nmols_ml:1*nmols_ml] = repeat([0.16, -0.16, -0.16, -0.16], outer =nx*ny)   # δr
    δr_ml[1 + 1*nmols_ml:2*nmols_ml] = repeat([0.0, 0.0, 0.0, 0.0], outer =nx*ny)
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

    flgs = zeros(Int64, 2*nmols_ml)
    flgs[1 + 0*nmols_ml:1*nmols_ml] = fill(1, nmols_ml)
    flgs[1 + 1*nmols_ml:2*nmols_ml] = fill(2, nmols_ml)

    # set the lower and upper limits
    lower, upper = low_high_limits(initial_state, flgs)

    # overlayer
    overlayer_states = zeros(Float64, 2*nmols_ol2)
    overlayer_states[1 + 0*nmols_ol2 : 1*nmols_ol2] = theta_ol[1:nmols_ol2]    # θ
    overlayer_states[1 + 1*nmols_ol2 : 2*nmols_ol2] = phi_ol[1:nmols_ol2]    # ϕ 

    # Declaration of a state vector including all dof to match other programs
    initial_state_all = zeros(Float64, ndofs_ml+4*nmols_ol2)
    initial_state_all[1 + ndofs_ml + 0*nmols_ol2 : ndofs_ml + 1*nmols_ol2] = theta_ol[1:nmols_ol2]    # θ
    initial_state_all[1 + ndofs_ml + 1*nmols_ol2 : ndofs_ml + 2*nmols_ol2] = phi_ol[1:nmols_ol2]    # ϕ 
    initial_state_all[1 + ndofs_ml + 2*nmols_ol2 : ndofs_ml + 4*nmols_ol2] = vec(δr_ol)   # δr

    initial_state_all[1 + 0*nmols_ml:2*nmols_ml] = initial_state
    initial_state_all[1 + 2*nmols_ml:5*nmols_ml] = vec(δr_ml)   # δr
    initial_state_all[ndofs_ml]                  = 10.0

    fig_ini = show_figure(initial_state_all, com0_ml, com0_ol, "Ininal", 0)
    # display(fig_ini)

####################
# Run optimization #
####################
# for z in -0.4:0.05:0.4
# Threads.@threads for i in 1:20
    # println(Threads.threadid())
    
    global δz_ml, δz_ol = 0.0, 10.0

    modified_state = random_coords(initial_state, flgs, [π, 2*π])

    # δr_ml[1 + 0*nmols_ml:1*nmols_ml] = zeros(nmols_ml)#+[0.0, 0.0, 0.0, 0.0]   # δr
    # δr_ml[1 + 1*nmols_ml:2*nmols_ml] = zeros(nmols_ml)#+[0.0, 0.0, 0.0, 0.0]
    δr_ml[1 + 2*nmols_ml:3*nmols_ml] = fill(δz_ml, nmols_ml)   # δr

    initial_state_all[1 + 0*nmols_ml:2*nmols_ml] = modified_state
    initial_state_all[1 + 2*nmols_ml:5*nmols_ml] = vec(δr_ml)   # δr
    initial_state_all[ndofs_ml]                  = δz_ol          # overlayer height deviation from c.-of-m.

    g_tol = 1e-8
    x_tol = 1e-8
    f_tol = 1e-8
    inner_optimizer = LBFGS()
    @time res = optimize(energy, lower, upper, modified_state, Fminbox(inner_optimizer), Optim.Options(g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, iterations = 2000))


################
# Show Results #
################
en_ini = energy(initial_state)
en_final = Optim.minimum(res)
# Set combined geometry
final_state = deepcopy(initial_state_all)
final_state[1 + 0*nmols_ml:2*nmols_ml] = res.minimizer     # θ
# write_to_file(joinpath(filepath, "3/x$z _$i.txt"), [[en_ini, en_final], [initial_state_all, final_state],[UInt8(Optim.converged(res))]]) 

# fig_ini = show_figure(initial_state_all, com0_ml, com0_ol, "Ininal", 0)
# fig_final = show_figure(final_state, com0_ml, com0_ol, "Final", 0)

# save(joinpath(filepath, "1/initial$i.png"), fig_ini)
# save(joinpath(filepath, "1/final$i.png"), fig_final)

# println("Initial energy: $en_ini")
# show_params(initial_state_all)
# println("Final energy: $en_final")
# show_params(final_state)
end
end