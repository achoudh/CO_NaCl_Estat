function test(values)
    val = copy(values)
    val[1] = 0
end

values = [1,2]
@time test(values)
println(values)