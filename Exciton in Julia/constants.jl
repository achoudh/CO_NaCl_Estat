ν0 = [2050.82, 2036.1]

μ00, μ11, μ01 = -0.112, -0.087, 0.105 # "Debyes"; μ00 and μ11: R.Disselkamp et al., Surface Science 240 (1990) 193-210; for 12C16O. μ01 calculated for 13C18O.
unit1 = 5034.12
unit2 = 7.51691023

# grid = parse(Int, readline("For N×N grid, N=? "))
#grid = 30
odown = 0.0
nx, ny = 8, 8
nmols_uc = 4
a0_surf = 3.99
nmols_ml = nmols_uc*nx*ny
grid = Int(sqrt(nmols_ml))

unitcell = [0.0 0.0 0.0;1.0 0.0 0.0;1.0 1.0 0.0;0.0 1.0 0.0]


θu = 37 * pi /180 #[42, 32, 32, 42] # 38 pm 3
θf = 124 * pi / 180 # 124 pm 4 # 22 pm 5
Φu = [12,58,12,58] .* pi/180# [0+Φu0, 90-Φu0, 90-Φu0, 0+Φu0]
Φf = [204, 261, 204, 261].* pi/180 # [180+Φu[1]-Φf0, 180+Φu[2]+Φf0, 180+Φu[3]+Φf0, 180+Φu[4]-Φf0]

#electric field
θe = 45.0*pi/180; # Tilt of incident beam
nar, ncr = 1.0, 1.52;
nrat = nar/ncr

Tx = 2*cos(θe)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2) / ((ncr/nar) + (cos(θe))^(-1)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2));
Ty = 2 / (1 + (ncr/nar)*(cos(θe))^(-1)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2));
Tz = 2*(sin(θe))^2 / (1 + (nar/ncr)*(cos(θe))^(-1)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2));

Tp, Ts = [[sqrt(Tx), 0.0, sqrt(Tz)],[0.0, sqrt(Ty), 0.0]]

ep, es = Tp, Ts;
