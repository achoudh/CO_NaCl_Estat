#include("../interactions.jl")

theta_m = fill(26.0,4)#[0.0,0.0,0.0,0.0]
phi_m = [0.0,0.0,180.0,180.0]
# theta_m = zeros(Float64,4)
# phi_m = zeros(Float64,4)
z_ml_ang = 3.1e-10
comml, bondlength_ml, ϕ_ml, θ_ml = monolayer(nx, ny, theta_m, phi_m)
comml1 = copy(comml)
z_ml = z_ml_ang/a0_surf

function opt_ml_full(x) #all coordinate of the xC_unitcell is set free except z_ml

    theta_ml = x[1+0*nmols_ml:1*nmols_ml] .* pi / 180.0 #fill(38.0,4) * pi / 180.0
    phi_ml = x[1+1*nmols_ml:2*nmols_ml] .* pi / 180
    ml_in = zeros(Float64, nmols_ml, 3)
    ml_in[:,1] = x[1+2*nmols_ml:3*nmols_ml]
    ml_in[:,2] = x[1+3*nmols_ml:4*nmols_ml]
    ml_in[:,3] = x[1+4*nmols_ml:5*nmols_ml] #fill(x[17],4)
    
    ml_in[:,1] = ml_in[:,1] - round.(ml_in[:,1])
    ml_in[:,2] = ml_in[:,2] - round.(ml_in[:,2])
    
    comml = comml1 + ml_in
    
    Vmlml = zeros(Float64,10)
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        local rvec12 = (comml[i,:] - comml[j,:])
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        Vmlml .+= site_site_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i],bondlength_ml[j], phi_ml[j], theta_ml[j])
    end
    
    Vmlsurf = zeros(Float64,3)
    for i in 1:nmols_ml
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = ml_in[i,:] .* a0_surf
        rvec[3] += z_ml_ang

        ml_o = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
        out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c, 4)
        Vmlsurf .+= [out_attr, out_rep,- out_disp]
    end

    out = Vmlml[1] + sum(Vmlsurf)

   return out*joule2wn/nmols_ml
end 

using Optim
x = zeros(Float64, 5*nmols_ml)
x[3] = 170
g_tol = 1e-6
x_tol = 1e-6
f_tol = 1e-6
res_mono = optimize(opt_ml_full, x, LBFGS(), Optim.Options(g_tol=g_tol,x_tol=x_tol,iterations=1000))
print(res_mono)
print(res_mono.minimizer)


#include("visualize.jl")
#display_structure_fullmono(res_mono.minimizer)






using Printf
function show_params(x)
        
    θ_ml = x[1+0*nmols_ml:1*nmols_ml]  #fill(38.0,4) * pi / 180.0
    ϕ_ml = x[1+1*nmols_ml:2*nmols_ml] 
    ml_in = zeros(Float64, nmols_ml, 3)
    ml_in[:,1] = x[1+2*nmols_ml:3*nmols_ml]
    ml_in[:,2] = x[1+3*nmols_ml:4*nmols_ml]
    ml_in[:,3] = x[1+4*nmols_ml:5*nmols_ml] #fill(x[17],4)

    out = hcat(θ_ml, ϕ_ml, ml_in[:,1].*a0_surf/1e-10, ml_in[:,2].*a0_surf/1e-10, (z_ml.+ml_in[:,3]).*a0_surf/1e-10)
    
    @printf("%-10s %-10s %-10s %-10s %-10s\n", "θ", "ϕ", "x", "y", "z")
    for i in 1:nmols_ml
    @printf("%f %f  %f  %f  %f\n",out[i,1],out[i,2],out[i,3],out[i,4],out[i,5])
    end
    return out
end


show_params(res_mono.minimizer)

using Plots
plotly()
function nacl()
    
    cls = 20
    nac = 10

    ll = nx*2 - 1
    mm = ny*2 - 1

    posna= [[i, j, 0] for i in 0:ll for j in 0:mm]
    poscl = [[i+0.5, j+0.5, 0] for i in 0:ll-1 for j in 0:mm-1]

    #na_x = [posna[i][1]  for i in eachindex(posna)]
    na_r = hcat(posna...)
    cl_r = hcat(poscl...)

    na = scatter3d(na_r[1,:], na_r[2,:], na_r[3,:], ms = nac, c=:violet)
    scatter3d!(cl_r[1,:], cl_r[2,:], cl_r[3,:], ms = cls, c=:green, alpha=0.8)
    return na
end


function display_structure_fullmono(x)
    theta_ml = x[1+0*nmols_ml:1*nmols_ml] .* pi / 180.0 #fill(38.0,4) * pi / 180.0
    phi_ml = x[1+1*nmols_ml:2*nmols_ml] .* pi / 180
    ml_in = zeros(Float64, nmols_ml, 3)
    ml_in[:,1] = x[1+2*nmols_ml:3*nmols_ml]
    ml_in[:,2] = x[1+3*nmols_ml:4*nmols_ml]
    ml_in[:,3] = x[1+4*nmols_ml:5*nmols_ml] #fill(x[17],4)
    vv, ww, bcc = [v, w, bc] ./a0_surf
    
    ml_o = zeros(Float64,nmols_ml,3)
    ml_c = zeros(Float64,nmols_ml,3)
    ml_bc = zeros(Float64,nmols_ml,3)
    comml = comml1 + ml_in
    
    for i in 1:nmols_ml
        rvec = comml[i,:]
        
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        ml_o[i,:] = rvec + [vv*stheta*cosphi, vv*stheta*sphi, vv*costheta]
        
        ml_c[i,:] = rvec + [-ww*stheta*cosphi, -ww*stheta*sphi, -ww*costheta]
        ml_bc[i,:] = rvec + [-bcc*stheta*cosphi, -bcc*stheta*sphi, -bcc*costheta]    
    end
    
    na = nacl()
    display(na)
    scatter3d!(ml_c[:,1],ml_c[:,2],ml_c[:,3],ms=9,c=:black,camera=(0,90,))
    scatter3d!(ml_o[:,1],ml_o[:,2],ml_o[:,3],ms=8,c=:red, ticks=nothing)
    
    return na
end

###








### Now start with a condition where the molecules are tilted.
function opt_ml_unit_full(x) #all coordinate of the xC_unitcell is set free except z_ml

    theta_m = x[1:4] .* pi / 180.0 #fill(38.0,4) * pi / 180.0
    phi_m = x[5:8] .* pi / 180
    ml_in = zeros(Float64, 4, 3)
    ml_in[:,1] = x[9:12]
    ml_in[:,2] = x[13:16]
    ml_in[:,3] = x[17:20]#fill(x[17],4)
    
    ml_in[:,1] = ml_in[:,1] - round.(ml_in[:,1])
    ml_in[:,2] = ml_in[:,2] - round.(ml_in[:,2])
    
    
    theta_ml = vec(repeat(theta_m[1:nmols_uc], outer=(1,nx*ny)))
    phi_ml = vec(repeat(phi_m[1:nmols_uc], outer=(1,nx*ny)))
    ml_in = repeat(ml_in, outer=(nx*ny))
    
    comml = comml1 + ml_in
    
    Vmlml = zeros(Float64,10)
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        local rvec12 = (comml[i,:] - comml[j,:])
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        Vmlml .+= site_site_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i],bondlength_ml[j], phi_ml[j], theta_ml[j])
    end
    
    Vmlsurf = zeros(Float64,3)
    for i in 1:nmols_ml
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = ml_in[i,:] .* a0_surf
        rvec[3] += z_ml_ang

        ml_o = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
        out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c, 4)
        Vmlsurf .+= [out_attr, out_rep,- out_disp]
    end

    out = Vmlml[1] + sum(Vmlsurf)

   return out*joule2wn/nmols_ml
end  


theta_m = [0.0,0.0,0.0,0.0]
phi_m = [0.0,0.0,0.0,0.0]
z_ml_ang = 3.1e-10
comml, bondlength_ml, ϕ_ml, θ_ml = monolayer(nx, ny, theta_m, phi_m)
comml1 = copy(comml)
x = zeros(Float64, 5*nmols_uc)
x[1:4] = theta_m
x[5:8] = phi_m
using Dates
st = now()
now()-st
g_tol = 1e-8
x_tol = 1e-8
f_tol = 1e-8
res_mono_2 = optimize(opt_ml_unit_full, x, LBFGS(), Optim.Options(g_tol=g_tol,x_tol=x_tol,f_tol=f_tol,iterations=1000))
print(res_mono_2)
print(res_mono_2.minimizer)


display_structure_fullmono(res_mono_2.minimizer)
na= display_structure_fullmono(x)
display(na)

y_loop = collect(3.2:0.1:5.5)
en_0 = []
for y in y_loop
    y = y*1e-10
    theta = [36.07, 37.44] .*pi/180
    phi = [-180, 180.0] .*pi/180
    com = [[0.0,0.0,0.0], [y/1.414, y/1.414, 0.0]]
    r12 = com[2] - com[1]
    vmlml = site_site_interaction(r12,v+w,phi[1],theta[1],v+w,phi[2],theta[2])
    push!(en_0, vmlml[1]*joule2wn)
end

plot(y_loop, en_0, c=:black, lw=3, label = "parallel", xlabel="distance/Å", ylabel="energy/cm-1")
plot!(y_loop, en_180,c=:black, ls=:dash, lw=3, label = "anti-parallel")
#default()

out_data = hcat(y_loop, en_0, en_180)

# out_file = open("twoCO.txt", "w")
# using DelimitedFiles
# println(out_file,"distance/Å\tparallel/cm-1\ta-parallel/cm-1")
# writedlm(out_file,out_data,"\t")
# close(out_file)