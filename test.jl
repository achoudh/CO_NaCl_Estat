###############################
# Now testing CO-CO potential #
###############################

using GLMakie

x = rand(5)
y = rand(5)
z = rand(5)
meshscatter(x,y,z, markersize = 0.1, color = :gray)



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