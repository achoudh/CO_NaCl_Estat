################################################################################
                                    ########
                                    # NOTE #
                                    ########

# This program finds the distance dependency of CO orientation
# Considers different PES for CO-NaCl


################################################################################

using Random
using LinearAlgebra
using Printf
using Optim
using GLMakie

#################
# Load programs #
#################

include("../constants.jl")
# include("lattice_construction.jl")
include("../visualization_makie.jl")
# include("ir_spectra.jl")
include("attraction_PES_arnab.jl")
include("../simulated_annealing.jl")

filepath = "C:/Users/achoudh/ownCloud/my work/CO_NaCl-estat/Estat_results"

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
            x_new = x[i] + δq[1]*r # As it should be 0 to 180
        elseif f == 2 # ϕ-type coordinates
            x_new = x[i] + δq[2]*r # As it should be 0 to 360
        elseif f == 3
            x_new = x[i] + (r - δq[3]) # As it should be -0.5 to 0.5
        elseif f == 4  # z coordinate
            x_new = x[i] + (r - δq[4])  # How much should it move
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

