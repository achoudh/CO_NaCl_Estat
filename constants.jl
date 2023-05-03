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
#z_ml = 3.35e-10/a0_surf # 2.6e-10   # monolayer-surface distance

v::Float64 = 0.4903e-10 #Oxygen
w::Float64 = 0.6437e-10 #Carbon
bc::Float64 = 0.0817e-10 # bond center
# Multipole constants
# see Jascha's Mathematica Notebook "P:\Research\CO-NaCl\Theory\Electrostatic model\MM_R_fit\R_MM_fit.nb"
rot::Matrix{Float64} = [cos(pi/4) sin(pi/4) 0.0; -sin(pi/4) cos(pi/4) 0.0; 0.0 0.0 1.0]

# dipole
a10::Float64 =  3.63272
a11::Float64 = -7.28409
a12::Float64 =  7.5581
a13::Float64 = -4.28056
a14::Float64 =  1.30636
a15::Float64 = -0.215764
a16::Float64 =  0.0182316
a17::Float64 = -0.000619188
# quadrupole
a20::Float64 = -14.5407
a21::Float64 =  28.2314
a22::Float64 = -25.6766
a23::Float64 =  11.9892
a24::Float64 =  -2.99969
a25::Float64 =   0.404491
a26::Float64 =  -0.0273732
a27::Float64 =   0.00071271
# octapole
a30::Float64 =  6.45176
a31::Float64 = -10.3863
a32::Float64 =  2.47699
a33::Float64 =  5.54757
a34::Float64 = -3.83351
a35::Float64 =  0.969093
a36::Float64 = -0.108525
a37::Float64 =  0.00450608

# Arrays of all odd integer combination (l,m) up to |l|,|m|<3
nlm::Int64 = 12
lodd::Vector{Int64} = [ 1, 1,-1,-1, 3, 3,-3,-3, 1, 1,-1,-1 ]
modd::Vector{Int64} = [ 1,-1,-1, 1, 1,-1,-1, 1, 3,-3,-3, 3 ]

global lm ::Vector{Float64} = [sqrt(lodd[i]^2 + modd[i]^2) for i in 1:nlm]
global lma::Vector{Float64} = [lm[i]*(-2*pi/a0_NaCl) for i in 1:nlm]
global ret_c::Vector{Float64} = [(eps_NaCl/a0_NaCl)*(1/lm[i])*(-1)^((lodd[i]+modd[i])/2) * (1/(1+exp(-pi*lm[i]))) for i in 1:nlm]

mom_C::Vector{Float64} = [0.18314*e, 0.33842*di_au, -0.90316*qu_au, -0.25179*oc_au, 0.13324*hx_au]
mom_O::Vector{Float64} = [-0.02320*e, -0.29304*di_au, 0.09203*qu_au, -0.09083*oc_au, -0.02669*hx_au]
mom_BC::Vector{Float64} = [-0.15994*e, 0.34962*di_au, 0.49745*qu_au, 0.45488*oc_au, 0.33552*hx_au]

# All constants below are reported as alpha0, alpha1, ro0,ro1, ro2
c_na_rep::Vector{Float64} = [4.5036e10, 0.4343e10, 2.9090e-10, -0.0636e-10, 0.0488e-10]
c_cl_rep::Vector{Float64} = [3.5542e10, 0.3156e10, 3.6047e-10, -0.0079e-10, 0.0948e-10]
o_na_rep::Vector{Float64} = [5.1882e10, -0.1221e10, 2.7192e-10, -0.0074e-10, -0.0455e-10]
o_cl_rep::Vector{Float64} = [3.8639e10, -0.0658e10, 3.2899e-10, 0.1018e-10, -0.0460e-10]
K_stone::Float64 = 4.3597482e-21
disp_coef::Vector{Float64} = [383.3, 256.6, 3935.9, 2633.0]/6.02214076*1e-80 #C-Na, O-Na, C-Cl, O-Cl

