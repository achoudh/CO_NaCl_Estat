

Rd = [[R0 * [xn[i] - xn[j] - nmols_ml * round((xn[i] - xn[j]) / nmols_ml),
              yn[i] - yn[j] - nmols_ml * round((yn[i] - yn[j]) / nmols_ml), 0]
       for j in 1:grid^2] for i in 1:grid^2] # Distance b/w any two molecules with peridic boundary conditions

rn = [[norm(Rd[i, j]) for j in 1:grid^2] for i in 1:grid^2] # |R|

en = [[ifelse(rn[i, j] == 0, [0, 0, 0], Rd[i, j] / rn[i, j])
        for j in 1:grid^2] for i in 1:grid^2] # unit vectors b/w any two molecules

# eu = [e_μs[pop[i, j] + 1] for i in 1:grid, j in 1:grid] # original code

