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
using CairoMakie

#################
# Load programs #
#################

include("../constants.jl")
# include("lattice_construction.jl")
include("../visualization_makie.jl")
# include("ir_spectra.jl")
# include("attraction_PES_arnab.jl")
# include("./site_surface_3DMA.jl")
include("./site_surface_arnab.jl")

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


function z_mono(x::Vector{Float64}, rvec::Vector{Float64})
    θ = x[1]*degrees
    ϕ = x[2]*degrees

    rvec[1] = rvec[1] - round(rvec[1])
    rvec[2] = rvec[2] - round(rvec[2])
    
    rvec *= a0_surf

    stheta, sphi, costheta, cosphi = sin(θ), sin(ϕ), cos(θ), cos(ϕ)
    unit_vec = [stheta*cosphi, stheta*sphi, costheta]

    ml_o  = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c  = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    
    out = site_surface_interaction(ml_in, unit_vec)[1] + mol_surf_rep_stone(ml_o, ml_c, 4)
    return out*joule2wn

end



zval = []
xval = collect(-5:0.1:5)
yval = collect(-5:0.1:5)
z = 3.35*1e-10/a0_surf


pot = zeros(Float64,size(xval,1),size(yval,1))
initial = zeros(Float64,size(xval,1),size(yval,1),2)
final = zeros(Float64,size(xval,1),size(yval,1),2)


ml_in = zeros(Float64,3)

for (i,x) in enumerate(xval)
for (j,y) in enumerate(xval)

    ml_in[1] = x
    ml_in[2] = y
    ml_in[3] = z

    xini = zeros(Float64,2)
    flgs = [1,2]

    initial_state = random_coords(xini, flgs, [π, 2π])
    low, high = low_high_limits(initial_state, flgs)
    
    initial[i,j,:] = initial_state

    g_tol = 1e-8
    x_tol = 1e-8
    f_tol = 1e-8
    
    inner_optimizer = LBFGS()
    res = optimize(x->z_mono(x, ml_in), low, high, initial_state, Fminbox(inner_optimizer), 
            Optim.Options(g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, iterations = 2000))
    
    pot[i,j] = minimum(res)
    final[i,j,:] = res.minimizer
    # push!(zval, z)
end
end 

#writedlm(joinpath(filepath, "singleCO_multiz_optim.txt"), hcat([zval, vcat(initial...)./degrees, vcat(final...)./degrees, pot]...))
display(heatmap(xval, yval, pot))


display(heatmap(xval, yval, final[:,:,1]))

final ./degrees

stheta, costheta = sin(final[:,:,1]), cos.(final[:,:,1])
sphi, cosphi = sin(final[:,:,2]), cos.(final[:,:,2])
unit_vec = [stheta.*cosphi, stheta.*sphi, costheta]

