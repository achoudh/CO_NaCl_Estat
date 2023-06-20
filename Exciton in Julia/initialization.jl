# Initialization

using Random
using PlotlyJS, CSV, DataFrames
using LinearAlgebra

#using Plots
include("functions.jl")
include("constants.jl")

com, θ, ϕ = monolayer(nx, ny, θu, Φu)
com2, θff, ϕff = monolayer(nx, ny, θf, Φf)

pop = zeros(Int,nmols_ml)

for i in eachindex(θ)
    if rand()<odown
        pop[i] = 1
    end
end
ev = fill(ν0[1], nmols_ml)
for i in eachindex(pop)
    if pop[i]==1
        θ[i], ϕ[i] = θff[i], ϕff[i]
        ev[i] = ν0[2]
    end
end

display_structure(com, θ, ϕ)

eu = zeros(Float64, nmols_ml, 3) #Orientation of the dipole moments of each vectors
for i in 1:nmols_ml
    eu[i,:] = [sin(θ[i]) * cos(ϕ[i]), sin(θ[i]) * sin(ϕ[i]), cos(θ[i])]
end


# define hstatT
hstatT = zeros(Float64, nmols_ml, nmols_ml)
for n1 = 1:nmols_ml, n2 = 1:nmols_ml
    if n1 == n2
        hstatT[n1,n2] = ev[n1] + unit1 * (μ11 - μ00) * μ00 * sum1(n1)
    else
        hstatT[n1,n2] = unit1 * μ01^2 * fij(n1, n2)
    end
end

# calculate eigenvalues and eigenvectors
eigenvals, eigenvecs = eigen(hstatT)

