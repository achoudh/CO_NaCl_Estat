include("..\\updated\\constants.jl")
include("..\\updated\\co_co_interactions.jl")
include("..\\updated\\site_surface_interaction.jl")
include("..\\updated\\lattice_construction.jl")
include("..\\updated\\co_co_Guo.jl")


comml, bondlength_ml, phi_ml, theta_ml = zeros(Float64, nx*ny, 3), zeros(Float64, nx*ny), zeros(Float64, nx*ny), zeros(Float64, nx*ny)
comol, bondlength_ol, phi_ol, theta_ol = zeros(Float64, nx*ny*nz, 3), zeros(Float64, nx*ny*nz), zeros(Float64, nx*ny*nz), zeros(Float64, nx*ny*nz)
rvec, ml_o, ml_c, ml_bc = zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3)


#global dz = 0.5
theta_m = zeros(Float64,4)
phi_m = zeros(Float64,4)
z_ml = 3.1e-10/a0_surf
comml, bondlength_ml, phi_ml, theta_ml = monolayer(theta_m, phi_m, z_ml)

Vmlml = 0.0 #zeros(Float64,1)
@time for i in 1:nmols_ml-1, j in i+1:nmols_ml
    local rvec12::Vector{Float64} = (comml[i,:] - comml[j,:])
    rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
    rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
    rvec12 = a0_surf .* rvec12
    # Vmlml .+= co_co_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i], bondlength_ml[j], phi_ml[j], theta_ml[j])
    Vmlml += co_co_int_nu(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
end
Vmlml



θ_over::Vector{Float64} = [pi*3.0/4.0, pi/4.0, 3.0*pi/4.0, pi/4.0, 3.0*pi/4.0, pi/4.0, 3.0*pi/4.0, pi/4.0]
ϕ_over::Vector{Float64} = [-pi/2.0, 0.0, -pi/2.0, 0.0, pi/2.0, pi, pi/2.0, pi]

comol, bondlength_ol, phi_ol, theta_ol = overlayer(θ_over, ϕ_over, dz, z_ml)
comml1 = copy(comml)
comol1 = copy(comol) 


function opt_ol_full(x::Vector{Float64}) #all coordinate of the xC_unitcell is set free except z_ml
    theta_m::Vector{Float64} = x[1:4] .* pi / 180.0 #fill(38.0,4) * pi / 180.0
    phi_m::Vector{Float64} = x[5:8] .* pi / 180
    ml_in = zeros(Float64, 4, 3)
    ml_in[:,1] = x[9:12]
    ml_in[:,2] = x[13:16]
    ml_in[:,3] = x[17:20]
    z_ml = 3.5e-10/a0_surf #+ x[21] #:5*nmols_ml] #fill(x[17],4)
    global dz = 0.5 #+ x[22] #x[2+4*nmols_ml]
    comml, bondlength_ml, phi_ml, theta_ml = monolayer(theta_m, phi_m, z_ml)
    
    ml_in[:,1] = ml_in[:,1] - round.(ml_in[:,1])
    ml_in[:,2] = ml_in[:,2] - round.(ml_in[:,2])

    comol, bondlength_ol, phi_ol, theta_ol = overlayer(θ_over, ϕ_over, z_ml, dz)
    comml1 = copy(comml)
    comol1 = copy(comol)

    ml_in = repeat(ml_in, outer=(nx*ny))
    comml = comml1 + ml_in
    
    #ol_in = zeros(Float64, nmols_ol, 3)
    Vmlml = 0.0 #zeros(Float64,1)
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        local rvec12 = (comml[i,:] - comml[j,:])
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        # Vmlml .+= co_co_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i], bondlength_ml[j], phi_ml[j], theta_ml[j])
        Vmlml += co_co_int_nu(rvec12, phi_ml[i], theta_ml[i], phi_ml[j], theta_ml[j])
    end
    
    Vmlsurf = zeros(Float64,3)
    for i in 1:nmols_ml
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = ml_in[i,:] .* a0_surf
        rvec[3] += comml1[i,3]*a0_surf

        ml_o = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
        out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c, 4)
        Vmlsurf .+= [out_attr, out_rep,- out_disp]
    end
    
    Vmlol = 0.0 #zeros(Float64, 11)
    for i in 1:nmols_ml, j in i:nmols_ol
        local rvec12 = (comml[i,:] .- comol[j,:])
        
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12

        # Vmlol .+= co_co_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i],bondlength_ol[j], phi_ol[j], theta_ol[j])
        
        Vmlol += co_co_int_nu(rvec12, phi_ml[i], theta_ml[i], phi_ol[j], theta_ol[j])
    end

    Volsurf = zeros(Float64,3)
    for i in 1:nmols_ol
        if (comol1[i,3] - dz*a0_CO/a0_surf -comml1[1,3]) <0.02
        rvec = [-0.5, -0.5, comol1[i,3]] .*a0_surf
        stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
        

        ml_o = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
        out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c, 4)
        Volsurf .+= [out_attr, out_rep,- out_disp]
        end
    end
    # out = (Vmlml[1]+Vmlol[1]+sum(Vmlsurf)+sum(Volsurf))*joule2wn/nmols_ml
    out = ((sum(Vmlsurf)+sum(Volsurf))*joule2wn + Vmlml + Vmlol)/nmols_ml
    
    return  out
end 



using Optim

x = zeros(Float64, 20)
@time opt_ol_full(x)

x[1:4] = [0.0,0.0,0.0,0.0]# fill(26.0,4)
g_tol = 1e-8

x_tol = 1e-8
f_tol = 1e-8
@time res_over = optimize(opt_ol_full, x, LBFGS(), Optim.Options(g_tol=g_tol,x_tol=x_tol, f_tol=f_tol, iterations = 2000))
print(res_over)
println(res_over.minimizer)
opt_ol_full(res_over.minimizer)

for i in 1:nmols_ol
    if abs(comol[i,3] - 0.5*a0_CO/a0_surf - comml[1,3]) <0.01
        println(comol[i,:])
    end
end

display_structure_unitmono(x)
display_structure_unitmono(res_over.minimizer)


# xfinal = res_over2.minimizer

# out_file = open("buried_1204.txt", "w")
# using DelimitedFiles
# println(out_file,minimum(res_over2))
# writedlm(out_file,xfinal,"\t")
# close(out_file)

minimum(res_over2)

#include("visualize.jl")
# p1 = display_structure_fullmono(res_mono.minimizer)

using Printf

function show_params(x)
    
    theta_m = x[1:4]  #fill(38.0,4) * pi / 180.0
    phi_m = x[5:8] 
    ml_in = zeros(Float64, 4, 3)
    ml_in[:,1] = x[9:12]
    ml_in[:,2] = x[13:16]
    #ml_in[:,3] = x[17:20]
    z_ml = fill(3.0e-10/a0_surf + x[21], 4)
    out = hcat(theta_m, phi_m, ml_in[:,1].*a0_surf/1e-10, ml_in[:,2].*a0_surf/1e-10, z_ml .*a0_surf/1e-10)
    
    @printf("%-10s %-10s %-10s %-10s %-10s\n", "θ", "ϕ", "x", "y", "z")
    for i in 1:4
    @printf("%f %f  %f  %f  %f\n",out[i,1],out[i,2],out[i,3],out[i,4],out[i,5])
    end
    return out
end

show_params(res_over2.minimizer)


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

    na = scatter3d(na_r[1,:], na_r[2,:], na_r[3,:], ms = nac, c=:violet, label = "Na")
    scatter3d!(cl_r[1,:], cl_r[2,:], cl_r[3,:], ms = cls, c=:green, alpha=0.8, label = "Cl")
    return na
end


function display_structure_unitmono(x)
    theta_m = x[1:4] .* pi / 180.0 #fill(38.0,4) * pi / 180.0
    phi_m = x[5:8] .* pi / 180
    ml_in = zeros(Float64, 4, 3)
    #ml_in[:,1] = x[9:12]
    #ml_in[:,2] = x[13:16]
    #ml_in[:,3] = x[17:20]
    global z_ml = 4.1e-10/a0_surf #+ x[17] #:5*nmols_ml] #fill(x[17],4)
    global dz = 0.5 #+ x[18] #x[2+4*nmols_ml]
    comml, bondlength_ml, phi_ml, theta_ml = monolayer(theta_m, phi_m, z_ml ) #fill(x[17],4)
    vv, ww, bcc = [v, w, bc] ./a0_surf
    
    ml_o = zeros(Float64,nmols_ml,3)
    ml_c = zeros(Float64,nmols_ml,3)
    ml_bc = zeros(Float64,nmols_ml,3)
    
    
    ml_in = repeat(ml_in, outer=(nx*ny))
    comml += ml_in
    
    for i in 1:nmols_ml
        rvec = comml[i,:]
        
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        ml_o[i,:] = rvec + [vv*stheta*cosphi, vv*stheta*sphi, vv*costheta]
        
        ml_c[i,:] = rvec + [-ww*stheta*cosphi, -ww*stheta*sphi, -ww*costheta]
        ml_bc[i,:] = rvec + [-bcc*stheta*cosphi, -bcc*stheta*sphi, -bcc*costheta]    
    end
    
    na = nacl()
    display(na)
    scatter3d!(ml_c[:,1],ml_c[:,2],ml_c[:,3],ms=9,c=:black,camera=(0,90,),label="C")
    scatter3d!(ml_o[:,1],ml_o[:,2],ml_o[:,3],ms=8,c=:red, ticks=nothing, label = "O")
    zlims!(0,10)
    return na
end


display_structure_unitmono(res_over.minimizer)


x_final = res_over2.minimizer
x = x_final
theta_m = x[1:4] .*pi/180 #fill(38.0,4) * pi / 180.0
phi_m = x[5:8] .*pi/180
ml_in = zeros(Float64, 4, 3)
ml_in[:,1] = x[9:12]
ml_in[:,2] = x[13:16]
ml_in[:,3] = x[17:20]
z_ml = 3.35e-10/a0_surf + x[21]#fill(3.0e-10/a0_surf + x[17], 4)

comml, bondlength_ml, phi_ml, theta_ml = monolayer(nx, ny, theta_m, phi_m, z_ml)
dz = 0.5+x[22]
comol, bondlength_ol, phi_ol, theta_ol = overlayer(nx, ny, nz, θ_over, ϕ_over, z_ml, dz)

ml_in = repeat(ml_in, outer=(nx*ny))
comml1 = copy(comml)
comol1 = copy(comol) 

# na, cl = nacl()
# c_ml, o_ml = display_CO(comml, θ_ml, ϕ_ml)

#display(plot([na, cl, c_ml, o_ml]))

using Random

pos = 3#round(Int,rand()*nmols_ml)

function flip(x)
    theta_ml[pos] = x[1]*pi/180
    phi_ml[pos] = x[2]*pi/180

    ml_in[pos,1:2] =x[3:4]
    ml_in[pos,3] = +0.06+x[5]
    ml_in[pos,1] = ml_in[pos,1] - round(ml_in[pos,1])
    ml_in[pos,2] = ml_in[pos,2] - round(ml_in[pos,2]) 
    comml = comml1 + ml_in
    
    #ol_in = zeros(Float64, nmols_ol, 3)
    Vmlml = zeros(Float64,10)
    for i in 1:nmols_ml-1, j in i+1:nmols_ml
        if i==pos
        local rvec12 = (comml[i,:] - comml[j,:])
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12
        Vmlml .+= site_site_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i],bondlength_ml[j], phi_ml[j], theta_ml[j])
        end
    end
    
    Vmlsurf = zeros(Float64,3)
    for i in 1:nmols_ml
        if i == pos
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        rvec = ml_in[i,:] .* a0_surf
        rvec[3] += comml1[i,3]*a0_surf

        ml_o = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
        out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c, 4)
        Vmlsurf .+= [out_attr, out_rep,- out_disp]
        end
    end
    
    Vmlol = zeros(Float64, 10)
    for i in 1:nmols_ml, j in i:nmols_ol
        if i==pos
        local rvec12 = (comml[i,:] .- comol[j,:])
        
        rvec12[1] = rvec12[1] - 2*nx*round(Int, rvec12[1]/(2*nx))
        rvec12[2] = rvec12[2] - 2*ny*round(Int, rvec12[2]/(2*ny))
        rvec12 = a0_surf .* rvec12

        Vmlol .+= site_site_interaction(rvec12, bondlength_ml[i], phi_ml[i], theta_ml[i],bondlength_ol[j], phi_ol[j], theta_ol[j])
        end
    end

    Volsurf = zeros(Float64,3)
    for i in 1:nmols_ol
        if i == pos
        if (comol1[i,3] - dz*a0_CO/a0_surf -comml1[1,3]) <0.02
        rvec = [-0.5, -0.5, comol1[i,3]] .*a0_surf
        stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
        

        ml_o = rvec + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
        ml_c = rvec + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
        ml_bc = rvec + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
        
        out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
        out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c, 4)
        Volsurf .+= [out_attr, out_rep,- out_disp]
        end
    end
    end
    out = (Vmlml[1]+Vmlol[1]+sum(Vmlsurf)+sum(Volsurf))
    return out*joule2wn/nmols_ml
end

x = zeros(Float64, 5)
x[1] = 150
g_tol, x_tol, f_tol = 1e-8, 1e-8, 1e-8
low = [-Inf, -Inf, -0.25, -0.25, -0.25]
high = [Inf, Inf, 0.25, 0.25, 0.25]

res_flip22 = optimize(flip, low, high, x, Fminbox(LBFGS()), Optim.Options(g_tol=g_tol,x_tol=x_tol,f_tol=f_tol, iterations = 1000))
# print(res_flip22)
print(res_flip22.minimizer)

x_flip = res_flip22.minimizer

theta_ml[pos] = x_flip[1]*pi/180
phi_ml[pos] = x_flip[2]*pi/180

ml_in[pos,1:3] =x_flip[3:5]


display_flip(comml, theta_ml, phi_ml, ml_in)


function display_flip(comml, theta_ml, phi_ml, ml_in)
    vv, ww, bcc = [v, w, bc] ./a0_surf
    ml_o = zeros(Float64,nmols_ml,3)
    ml_c = zeros(Float64,nmols_ml,3)
    ml_bc = zeros(Float64,nmols_ml,3)

    comml = comml + ml_in

    for i in 1:nmols_ml
        rvec = comml[i,:]
        
        stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
        ml_o[i,:] = rvec + [vv*stheta*cosphi, vv*stheta*sphi, vv*costheta]
        
        ml_c[i,:] = rvec + [-ww*stheta*cosphi, -ww*stheta*sphi, -ww*costheta]
        ml_bc[i,:] = rvec + [-bcc*stheta*cosphi, -bcc*stheta*sphi, -bcc*costheta]    
    end

    na = nacl()
    display(na)
    scatter3d!(ml_c[:,1],ml_c[:,2],ml_c[:,3],ms=9,c=:black,camera=(0,90,),label="C")
    scatter3d!(ml_o[:,1],ml_o[:,2],ml_o[:,3],ms=8,c=:red, ticks=nothing, label = "O")
    return na
end