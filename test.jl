###############################
# Now testing CO-CO potential #
###############################

include("visualization_makie.jl")

data = readdlm("buried_ov_fixed_dof.txt")

no_of_iteration = data[2] 
energies = [i for i in data[4,:] if i != ""]
dof = data[6]
states = [[i for i in data[j,:] if i != ""] for j in 7:7+no_of_iteration]

accepted_step_energies = [i for i in data[7+no_of_iteration+3,:] if i != ""] 


fig = Figure(resolution = (600, 600))
ax1 = set_axis(fig[1,1], L"Energie/cm$^{-1}$", "")
hist!(ax1, energies[2:end], 10)
ax2 = set_axis(fig[2,1], "Steps", L"Energie/cm$^{-1}$")
lines!(ax2, accepted_step_energies)
display(fig)

display(GLMakie.Screen(),show_figure(states[1], com0_ml, com0_ol, "Initial"))
for (i,s) in enumerate(states[2:end])
    display(GLMakie.Screen(),show_figure(s, com0_ml, com0_ol, ""))
end



###############################
# Induction #
###############################



# CO struct with multipoles and polarizabilities
struct CO
    q::Vector{Float64} # multipole moments 
    α::Matrix{Float64} # polarizability tensor
end

# Compute induced dipole at site i on mol1 due to fields from mol2 
function induce(mol1::CO, mol2::CO, i::Int)
    Ei = electric_field(mol2, mol1.q[i]) # Ewald field at i from mol2
    μi = mol1.α[i,:] * Ei # induced dipole
    return μi 
end

# Electric field at site i from periodic multipoles 
function electric_field(mol::CO, i::Int)
    Ei = zeros(3)
    for j in eachindex(mol.q)
        r = mol[j] - mol[i]
        Ei += efield(mol.q[j], r) # from Ewald sum
    end
    return Ei
end 

# Induction energy between two CO molecules
function induction(mol1::CO, mol2::CO)
    E = 0.0
    for i in eachindex(mol1.q)
        μi = induce(mol1, mol2, i) # induced dipole on mol1
        E -= 0.5 * μi' * electric_field(mol2, μi) # interaction energy
    end
    return E
end