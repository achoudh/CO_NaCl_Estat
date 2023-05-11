# Physical and mathematical constants
#pi = acos(-1.0)
sqrt2::Float64 = sqrt(2.0)
eps0::Float64 = 8.854187817e-12            # vacuum permittivity ( C/(Vm) = F/m = A^2 s^4 / (kg m^3) )
eps4pi::Float64 = 1.0/(4.0*pi*eps0)
bohr_radius::Float64 = 5.2917721067e-11    # m
boltzmann::Float64 = 1.380649e-23
di_au::Float64 = 8.478353551956304e-30     # e*a0   (C m)
qu_au::Float64 = 4.486551483698324e-40     # e*a0^2 (C m^2)
oc_au::Float64 = 2.374180799670829e-50     # e*a0^3 (C m^3)
hx_au::Float64 = oc_au*bohr_radius
nu::Int64 = 1                          # valency of NaCl
e::Float64 = 1.6021766208e-19              # elementary charge (C)
eps_NaCl::Float64 = nu*e/(pi*eps0)


# Unit conversion factors
joule2wn::Float64 = 1.0/1.986445824171758e-23
m2au::Float64 = 1.0/bohr_radius
degrees::Float64 = pi/180.0

# CO layer parameters
nmols_uc::Int64 = 4
nmols_ucol::Int64 = 8
nx::Int64 = 2
ny::Int64 = 2
nz::Int64 = 2

nmols_ml::Int64 = nmols_uc*nx*ny
nmols_ol::Int64 = nmols_ucol*nx*ny*nz


r_CO::Float64 = 1.14e-10       # CO bondlength in m
a0_CO::Float64 = 5.64e-10      # CO layer lattice constant, m
a0_NaCl::Float64 = 5.64e-10    # NaCl lattice constant, m
a0_surf::Float64 = 3.99e-10    # NaCl surface lattice constant
# NaCl lattice vectors
b1::Vector{Float64} = [0,1,0]*a0_surf
b2::Vector{Float64} = [1,0,0]*a0_surf
b3::Vector{Float64} = [0,0,1]*a0_NaCl
# Na and Cl positions for the 1st and 2nd layers
nl_surf = 2
pos_surf::Vector{Matrix{Float64}} = [ [[0.0, 0.0, 0.0] [-0.5*a0_surf, -0.5*a0_surf, 0.0]],
             [[-0.5*a0_surf, -0.5*a0_surf, 0.0]-b3 [0.0, 0.0, 0.0]-b3] ]
#[push!(pos_surf, reverse(pos_surf[end] .- b3,dims=2) ) for l in 1:nl_surf]
#z_ml = 3.35e-10/a0_surf # 2.6e-10   # monolayer-surface distance

v::Float64 = 0.4903e-10 #Oxygen
w::Float64 = 0.6437e-10 #Carbon
bc::Float64 = 0.0817e-10 # bond center

rot::Matrix{Float64} = [cos(pi/4) sin(pi/4) 0.0; -sin(pi/4) cos(pi/4) 0.0; 0.0 0.0 1.0]

# CO-CO NN-PES data
# Dimensions of arrays for NN-PES
s0_nn = 7
s1_nn = 45
s2_nn = 45

w1_nn::Matrix{Float64} = transpose(readdlm("w1.txt"))
w2_nn::Matrix{Float64} = transpose(readdlm("w2.txt"))
b1_nn::Matrix{Float64} = readdlm("b1.txt")
b2_nn::Matrix{Float64} = readdlm("b2.txt")
w3_nn::Matrix{Float64} = readdlm("w3.txt")
rg_nn::Matrix{Float64} = transpose(readdlm("rg.txt"))
vg_nn::Vector{Float64} = [-702.59000000000003, 22446.305300000000]
b3_nn::Float64 = 2.6335296521641700



# Arrays of all odd integer combination (l,m) up to |l|,|m|<3
nlm::Int64 = 12
lodd::Vector{Int64} = [ 1, 1,-1,-1, 3, 3,-3,-3, 1, 1,-1,-1 ]
modd::Vector{Int64} = [ 1,-1,-1, 1, 1,-1,-1, 1, 3,-3,-3, 3 ]

lm ::Vector{Float64} = sqrt.(lodd.^2 + modd.^2)
lma::Vector{Float64} = -2*pi/a0_NaCl .* lm
ret_c::Vector{Float64} = eps_NaCl/a0_NaCl ./ lm .* (-1).^((lodd+modd)/2) ./ (1 .+ exp.(-pi*lm))

mom_C::Vector{Float64} = [0.18314*e, 0.33842*di_au, -0.90316*qu_au, -0.25179*oc_au, 0.13324*hx_au]
mom_O::Vector{Float64} = [-0.02320*e, -0.29304*di_au, 0.09203*qu_au, -0.09083*oc_au, -0.02669*hx_au]
mom_BC::Vector{Float64} = [-0.15994*e, 0.34962*di_au, 0.49745*qu_au, 0.45488*oc_au, 0.33552*hx_au]

# All constants below are reported as alpha0, alpha1, ro0,ro1, ro2
c_na_rep::Vector{Float64} = [4.5036e10, 0.4343e10, 2.9090e-10, -0.0636e-10, 0.0488e-10]
c_cl_rep::Vector{Float64} = [3.5542e10, 0.3156e10, 3.6047e-10, -0.0079e-10, 0.0948e-10]
o_na_rep::Vector{Float64} = [5.1882e10, -0.1221e10, 2.7192e-10, -0.0074e-10, -0.0455e-10]
o_cl_rep::Vector{Float64} = [3.8639e10, -0.0658e10, 3.2899e-10, 0.1018e-10, -0.0460e-10]
rep_coeffs = [ [c_na_rep, o_na_rep], [c_cl_rep, o_cl_rep] ]
K_stone::Float64 = 4.3597482e-21

# Dispersion coefficients
# [ [C-Na, O-Na], [C-Cl, O-Cl] ]
disp_coef::Matrix{Float64} = [ [383.3 256.6]; [3935.9 2633.0] ]/6.02214076*1e-80 



