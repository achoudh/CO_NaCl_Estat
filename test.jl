function test(values)
    val = copy(values)
    val[1] = 0
end

values = [1,2]
@time test(values)


println(values)


println("mol\tθ°\t", "ϕ°" )
for i in 1:8
println(i,"\t",θ_uc[i]/degrees,"\t", ϕ_uc[i]/degrees)
end