using LinearAlgebra

include("constants.jl")


function fit_CO_CO_all(zfit)
    include("co_co_Guo.jl")
    eml = []
    eml2 = []
    eml3 = []
    eml4 = []
    eml5 = []
    for y in zfit
        y = y*1e-10
        theta1 = [0.0, 180.0] *pi/180 
        theta2 = [0.0,   0.0] *pi/180
        phi = [0.0, 0.0] *pi/180 
        theta3 = [37.0, 37.0] *pi/180 
        phi3 = [-90, 90] *pi/180
        phi4 = [90, 90] *pi/180 
        theta5 = [37.0, 143.0] *pi/180
        phi5 = [-90, -90] *pi/180
        
        com = [[0.0, 0.0, 0.0], [0.0,  y, 0.0]]
        r12 = com[1] - com[2]
        vmlml1 = co_co_interaction(r12, phi[1],theta1[1],phi[2],theta1[2])
        vmlml2 = co_co_interaction(r12, phi[1],theta2[1],phi[2],theta2[2])
        vmlml3 = co_co_interaction(r12,phi3[1],theta3[1],phi3[2],theta3[2])
        vmlml4 = co_co_interaction(r12,phi4[1],theta3[1],phi4[2],theta3[2])
        vmlml5 = co_co_interaction(r12,phi5[1],theta5[1],phi5[2],theta5[2])
        

        push!(eml , vmlml1[1])
        push!(eml2, vmlml2[1])
        push!(eml3, vmlml3[1])
        push!(eml4, vmlml4[1])
        push!(eml5, vmlml5[1])
    end
    return eml, eml2, eml3, eml4, eml5
end 


function fit_CO_CO_all_2(zfit)
    include("co_co_interactions.jl")
    eml = []
    eml2 = []
    eml3 = []
    eml4 = []
    eml5 = []
    for y in zfit
        y = y*1e-10
        theta1 = [0.0, 180.0] *pi/180 
        theta2 = [0.0,   0.0] *pi/180
        phi = [0.0, 0.0] *pi/180 
        theta3 = [37.0, 37.0] *pi/180 
        phi3 = [-90, 90] *pi/180
        phi4 = [90, 90] *pi/180 
        theta5 = [37.0, 143.0] *pi/180
        phi5 = [-90, -90] *pi/180
        
        com = [[0.0, 0.0, 0.0], [0.0,  y, 0.0]]
        r12 = com[1] - com[2]
        vmlml1 = co_co_interaction(r12, phi[1],theta1[1],phi[2],theta1[2])
        vmlml2 = co_co_interaction(r12, phi[1],theta2[1],phi[2],theta2[2])
        vmlml3 = co_co_interaction(r12,phi3[1],theta3[1],phi3[2],theta3[2])
        vmlml4 = co_co_interaction(r12,phi4[1],theta3[1],phi4[2],theta3[2])
        vmlml5 = co_co_interaction(r12,phi5[1],theta5[1],phi5[2],theta5[2])
        

        push!(eml , vmlml1[1])
        push!(eml2, vmlml2[1])
        push!(eml3, vmlml3[1])
        push!(eml4, vmlml4[1])
        push!(eml5, vmlml5[1])
    end
    return eml, eml2, eml3, eml4, eml5
end 

zval = collect(2.5:4.5)
println(fit_CO_CO_all(zval))
println(fit_CO_CO_all_2(zval))
