###############################
# Now testing CO-CO potential #
###############################

using GLMakie

x = rand(5)
y = rand(5)
z = rand(5)
meshscatter(x,y,z, markersize = 0.1, color = :gray)
