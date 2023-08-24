# CO_NaCl_Estat
Simulated annealing Metropolis search for the global minimum 
for CO monolayer at NaCl(100)

    CO-NaCl interaction is modeled by Meredith & Stone potential

    CO-CO interaction is modeled by NN-PES by Guo's group

    IR Spectra is generated using Exciton model based on R. Disselkamp, H. Chang, and G.E. Ewing, Surf. Sci. 240 (1990) 193-210



####################################
Current update            24-08-2023      

1. Steps and initial conditions are 
2. limited to only have the C-down  


in line 16 of simulated annealing, a new line is added to make sure that CO is always C-down
         x_new = x_new > pi/2 ?  pi - x_new : x_new  

####################################

