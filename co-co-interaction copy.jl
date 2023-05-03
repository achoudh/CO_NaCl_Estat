using CSV, DataFrames
using DelimitedFiles


include("constants.jl")
include("co_co_interactions.jl")
include("../Plots_default.jl")

#%%
### Now my own interactions
# define the construction
θ = [134.95, 45.05] *pi/180
ϕ = [-90, 90] *pi/180
θ = [65.2,114.8] *pi/180
ϕ = [90, -90] *pi/180
com = [[0.0 0.0 0.0];[0.0 2.6 0.0]]

ml_o = zeros(Float64,2,3)
ml_c = zeros(Float64,2,3)
vv = v*1e10
ww = w*1e10
for i in 1:2
stheta, sphi, costheta, cosphi = sin(θ[i]), sin(ϕ[i]), cos(θ[i]), cos(ϕ[i])
ml_o[i,:] .= com[i,:] .+ [vv*stheta*cosphi, vv*stheta*sphi, vv*costheta]
ml_c[i,:] .= com[i,:] .+ [-ww*stheta*cosphi, -ww*stheta*sphi, -ww*costheta]
end

scatter3d(ml_o[:,1], ml_o[:,2], ml_o[:,3], c=:red, ms=8, xlabel="x", ylabel="y", zlabel="z")
scatter3d!(ml_c[:,1], ml_c[:,2], ml_c[:,3], c=:black, ms=9, ticks = nothing)
scatter3d!(com[:,1], com[:,2], com[:,3], c=:blue, ms=3, ticks = nothing)
#scatter3d!(ml_c[:,1], ml_c[:,2], ml_c[:,3], c=:black, ms=9, ticks = nothing)

xlims!(-1,1)
zlims!(-4,4)

lim1 = 3.2
lim2 = 5.0

### Import data and plot
data_CO_CO_2 = CSV.read("..\\New_analysis\\CO-CO\\ener_ccsdt_config2.dat", DataFrame, delim=" ", header=2) 
z_CO_2 = data_CO_CO_2[:,1]
en_CO_2 = data_CO_CO_2[:,4]
z_fit_2 = z_CO_2[argmin([abs(z - lim1) for z in z_CO_2]): argmin([abs(z - lim2) for z in z_CO_2])]
en_fit_2 = en_CO_2[argmin([abs(z - lim1) for z in z_CO_2]): argmin([abs(z - lim2) for z in z_CO_2])]


#%% for config 1 
com1 = [7.812825, 1.284041, 18.278374]
com2 = [7.879290, 0.614623, 19.201797]
comc2 = [3.802138, 2.925753, 18.265551]

vec = comc2 - com1
mag = norm(vec)
vec1 = vec / mag

data_CO_CO_1 = CSV.read("..\\New_analysis\\CO-CO\\energies_config1.dat", DataFrame, delim=" ", header=2) 
z_CO_1 = data_CO_CO_1[:,1]
en_CO_1 = data_CO_CO_1[:,4]


data_CO_CO_par = CSV.read("..\\New_analysis\\CO-CO\\ener_parallel_straight.dat", DataFrame, delim=" ", header=2) 
z_CO_par = data_CO_CO_par[:,1]
en_CO_par = data_CO_CO_par[:,4]
z_fit_par = z_CO_par[argmin([abs(z - lim1) for z in z_CO_par]): argmin([abs(z - lim2) for z in z_CO_par])]
en_fit_par = en_CO_par[argmin([abs(z - lim1) for z in z_CO_par]): argmin([abs(z - lim2) for z in z_CO_par])]


data_CO_CO_anti_par = CSV.read("..\\New_analysis\\CO-CO\\ener_anti_parallel_straight.dat", DataFrame, delim=" ", header=2) 
z_CO_anti_par = data_CO_CO_anti_par[:,1]
en_CO_anti_par = data_CO_CO_anti_par[:,4]
z_fit_anti_par = z_CO_anti_par[argmin([abs(z - lim1) for z in z_CO_anti_par]): argmin([abs(z - lim2) for z in z_CO_anti_par])]
en_fit_anti_par = en_CO_anti_par[argmin([abs(z - lim1) for z in z_CO_anti_par]): argmin([abs(z - lim2) for z in z_CO_anti_par])]


data_CO_CO_anti_par_tilt = CSV.read("..\\New_analysis\\CO-CO\\ener_antiparallel.dat", DataFrame, delim=" ", header=2) 
z_CO_anti_par_tilt = data_CO_CO_anti_par_tilt[:,1]
en_CO_anti_par_tilt = data_CO_CO_anti_par_tilt[:,4]
z_fit_anti_par_tilt = z_CO_anti_par_tilt[argmin([abs(z - lim1) for z in z_CO_anti_par_tilt]): argmin([abs(z - lim2) for z in z_CO_anti_par_tilt])]
en_fit_anti_par_tilt = en_CO_anti_par_tilt[argmin([abs(z - lim1) for z in z_CO_anti_par_tilt]): argmin([abs(z - lim2) for z in z_CO_anti_par_tilt])]



data_CO_CO_par_tilt = CSV.read("..\\New_analysis\\CO-CO\\ener_parallel_tilt.dat", DataFrame, delim=" ", header=2) 
z_CO_par_tilt = data_CO_CO_par_tilt[:,1]
en_CO_par_tilt = data_CO_CO_par_tilt[:,4]
z_fit_par_tilt = z_CO_par_tilt[argmin([abs(z - lim1) for z in z_CO_par_tilt]): argmin([abs(z - lim2) for z in z_CO_par_tilt])]
en_fit_par_tilt = en_CO_par_tilt[argmin([abs(z - lim1) for z in z_CO_par_tilt]): argmin([abs(z - lim2) for z in z_CO_par_tilt])]
###


function fit_CO_CO_all(k)
    global k21, k23 = 1.0, 1.0 # k21, k22, k23
    #k11, k12, k13 = #fill(1.0,3)
    #global k21, k22, k23 = fill(0.65,3)
    #lobal k11, k12, k13 = [k11, k12, k13] .*1e10
    
    global k11 = 1.0
    global k13 = 1.0

    eml = []
    eml2 = []
    eml3 = []
    eml4 = []
    eml5 = []
    for y in z_fit_par
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
        vmlml1 = co_co_interaction(r12,v+w,phi[1],theta1[1],v+w,phi[2],theta1[2])
        vmlml2 = co_co_interaction(r12,v+w,phi[1],theta2[1],v+w,phi[2],theta2[2])
        vmlml3 = co_co_interaction(r12,v+w,phi3[1],theta3[1],v+w,phi3[2],theta3[2])
        vmlml4 = co_co_interaction(r12,v+w,phi4[1],theta3[1],v+w,phi4[2],theta3[2])
        vmlml5 = co_co_interaction(r12,v+w,phi5[1],theta5[1],v+w,phi5[2],theta5[2])
        

        push!(eml , vmlml1[1]*joule2wn)
        push!(eml2, vmlml2[1]*joule2wn)
        push!(eml3, vmlml3[1]*joule2wn)
        push!(eml4, vmlml4[1]*joule2wn)
        push!(eml5, vmlml5[1]*joule2wn)
    end
    res1 = sqrt(sum((en_fit_anti_par - eml) .^2))
    res2 = sqrt(sum((en_fit_par - eml2) .^2))
    res3 = sqrt(sum((en_fit_2 - eml3) .^2))
    res4 = sqrt(sum((en_fit_par_tilt - eml4) .^2))
    res5 = sqrt(sum((en_fit_anti_par_tilt - eml5) .^2))
    return res5 + res1 + res2#+ 
end 


function fit_CO_CO_all_2(zfit)
    
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
        vmlml1 = co_co_interaction(r12,v+w,phi[1],theta1[1],v+w,phi[2],theta1[2])
        vmlml2 = co_co_interaction(r12,v+w,phi[1],theta2[1],v+w,phi[2],theta2[2])
        vmlml3 = co_co_interaction(r12,v+w,phi3[1],theta3[1],v+w,phi3[2],theta3[2])
        vmlml4 = co_co_interaction(r12,v+w,phi4[1],theta3[1],v+w,phi4[2],theta3[2])
        vmlml5 = co_co_interaction(r12,v+w,phi5[1],theta5[1],v+w,phi5[2],theta5[2])
        

        push!(eml , vmlml1[1]*joule2wn)
        push!(eml2, vmlml2[1]*joule2wn)
        push!(eml3, vmlml3[1]*joule2wn)
        push!(eml4, vmlml4[1]*joule2wn)
        push!(eml5, vmlml5[1]*joule2wn)
    end
    return eml, eml2, eml3, eml4, eml5
end 


function fit_CO_CO(zfit)

    en_0 = []
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
        vmlml1 = co_co_interaction(r12,v+w,phi[1],theta1[1],v+w,phi[2],theta1[2])
        vmlml2 = co_co_interaction(r12,v+w,phi[1],theta2[1],v+w,phi[2],theta2[2])
        vmlml3 = co_co_interaction(r12,v+w,phi3[1],theta3[1],v+w,phi3[2],theta3[2])
        vmlml4 = co_co_interaction(r12,v+w,phi4[1],theta3[1],v+w,phi4[2],theta3[2])
        vmlml5 = co_co_interaction(r12,v+w,phi5[1],theta5[1],v+w,phi5[2],theta5[2])
        push!(en_0, vmlml1[1]*joule2wn)
    end
    return en_0    
end


function fit_CO_CO_nu(zfit)

    #global k21, k22, k23 = fill(0.65,3)
    #global k1 = 7
    #global k11, k12, k13 = [k11, k12, k13] .*1e10

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
        #theta3 = [45.4, 45.4] *pi/180 
        
        phi3 = [90, -90] *pi/180
        phi4 = [90, 90] *pi/180 
        theta5 = [37.0, 143.0] *pi/180
        phi5 = [-90, -90] *pi/180
        
        com = [[0.0, 0.0, 0.0], [0.0,  y, 0.0]]
        r12 = com[1] - com[2]
        vmlml1 = co_co_int_nu(r12,phi[1],theta1[1],phi[2],theta1[2])
        vmlml2 = co_co_int_nu(r12,phi[1],theta2[1],phi[2],theta2[2])
        vmlml3 = co_co_int_nu(r12,phi3[1],theta3[1],phi3[2],theta3[2])
        vmlml4 = co_co_int_nu(r12,phi4[1],theta3[1],phi4[2],theta3[2])
        vmlml5 = co_co_int_nu(r12,phi5[1],theta5[1],phi5[2],theta5[2])
        

        push!(eml , vmlml1[1])
        push!(eml2, vmlml2[1])
        push!(eml3, vmlml3[1])
        push!(eml4, vmlml4[1])
        push!(eml5, vmlml5[1])
    end
    return eml, eml2, eml3, eml4, eml5
end



z_data = hcat(z_fit_anti_par, z_fit_par, z_fit_2, z_fit_par_tilt, z_fit_anti_par_tilt)
en_data = hcat(en_fit_anti_par, en_fit_par, en_fit_2, en_fit_par_tilt, en_fit_anti_par_tilt)



theta1 = [0.0, 180.0] *pi/180
phi1 = [0.0, 0.0] *pi/180

com = [[0.0, 0.0, 0.0], [0.0, 3.6181e-10,  0.0]]
r12 = com[1] - com[2]
vmlml1 = co_co_interaction(r12,v+w,phi1[1],theta1[1],v+w,phi1[2],theta1[2])



p4 = plot(z_CO_par_tilt, en_CO_par_tilt)
p4 = scatter!(z_CO_par_tilt,fit_CO_CO_all_2(z_CO_par_tilt)[4])
p4 = plot!(z_CO_par_tilt,fit_CO_CO_nu(z_CO_par_tilt)[4], color=:black)
display(p4)

p2 = scatter(z_CO_par, en_CO_par)
p2 = plot!(z_CO_par,fit_CO_CO_all_2(z_CO_par)[2])
p2 = plot!(z_CO_par,fit_CO_CO_nu(z_CO_par)[2], color=:black)
display(p2)

p1 = scatter(z_CO_anti_par, en_CO_anti_par)
p1 = plot!(z_CO_anti_par,fit_CO_CO_all_2(z_CO_anti_par)[1])
plot!(z_CO_anti_par,fit_CO_CO_nu(z_CO_anti_par)[1], color=:black)
display(p1)

p3 = scatter(z_CO_2, en_CO_2)
p3 = plot!(z_CO_2,fit_CO_CO_all_2(z_CO_2)[3])
p3 = plot!(z_CO_2,fit_CO_CO_nu(z_CO_2)[3], color = :black)
display(p3)

p5 = scatter(z_CO_anti_par_tilt, en_CO_anti_par_tilt)
p5 = plot!(z_CO_anti_par_tilt,fit_CO_CO_all_2(z_CO_anti_par_tilt)[5])
p5 = plot!(z_CO_anti_par_tilt,fit_CO_CO_nu(z_CO_anti_par_tilt)[5], color = :black)
display(p5)

include("co_co_Guo.jl")

fit_CO_CO_nu(z_CO_anti_par_tilt)[1]

plotly()


theta2 = [0.0,   0.0] *pi/180
phi = [0.0, 0.0] *pi/180 

com = [[0.0, 0.0, 0.0], [0.0, 3.6181e-10,  0.0]] #4.3335e-10,
r12 = com[1]-com[2]

vmlml2 = co_co_interaction(r12,v+w,phi[1],theta2[1],v+w,phi[2],theta2[2])[1] *joule2wn