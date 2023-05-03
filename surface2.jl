
#%% create surface below molecule 
# b1 = [0,1,0]
# b2 = [1,0,0]

# na = []
# cl = []

# k = 2
# for i in -k:k, j in -k:k
#     push!(na, [0.0, 0.0,0.0] + i .*b1 + j.*b2)
#     push!(cl, [-0.5,-0.5,0.0] + i .*b1 + j.*b2)
# end

# na = hcat(na...)'
# cl = hcat(cl...)'

# using PlotlyJS

# nap = scatter3d(x=na[:,1], y=na[:,2], z=zeros(Float64,length(na)), mode="markers", marker=attr(size=10, color="green", opacity = 0.8))
# clp = scatter3d(x=cl[:,1], y=cl[:,2], z=zeros(Float64,length(na)), mode="markers", marker_size=20, marker_color ="blue")

# p = plot([nap, clp])
# display(p)

#%% See the impact of size of surface on repulsion
# include("interactions.jl")
# include("constants.jl")

# ml_in = [0.1, 0.1, 2.6*1e-10/a0_surf]

# ml_in *= a0_surf

# en =[]

# for k in 0:6
# ml_in = [0.1, 0.1, 2.6*1e-10/a0_surf]
# ml_in *= a0_surf

# push!(en, mol_surf_rep(ml_in,0,0,0,k)*joule2wn)

# end
# pic = scatter(x=collect(0:6),y= en, marker=(color ="black", size = 10))
# layout = Layout(xaxis_title="size of surface", yaxis_title="Energy/cm-1", title="Z = 2.6 Ang")
# plot(pic,layout)
include("constants.jl")
include("interactions.jl")
include("New_analysis/site_surface_interaction.jl")
#using Printf
using Optim

function z_mono2(x)
    theta = x[1] *pi/180
    #theta1 = 32 *pi/180
    phi = x[2] *pi/180
    ml_in = zeros(Float64,3)
    
    ml_in[1:2] = x[3:4]

    ml_in[1] = ml_in[1] - round(ml_in[1])
    ml_in[2] = ml_in[2] - round(ml_in[2])
    ml_in .*= a0_surf
    ml_in[3] = (3.1 + x[5])*1e-10
    stheta, sphi, costheta, cosphi = sin(theta), cos(phi), cos(theta), cos(phi)
    
    ml_o = ml_in + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c = ml_in + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ml_bc = ml_in + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
    out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c,4)

    return (out_attr+out_rep-out_disp)*joule2wn
end

x = zeros(Float64,5)
x[1] = 180
res2 = optimize(z_mono2, x, LBFGS(), Optim.Options(x_tol = 1e-8,f_tol = 1e-8,iterations = 50))
print(res2)
print(res2.minimizer)
z_mono2(x)

en = []
for i in 0:0.1:1.5
    theta = 0.0*pi/180
    #theta1 = 32 *pi/180
    phi = 0.0   #x[2] *pi/180
    ml_in = zeros(Float64,3)
    
    #ml_in[1:2] = x[3:4]

    #ml_in[1] = ml_in[1] - round(ml_in[1])
    #ml_in[2] = ml_in[2] - round(ml_in[2])
    ml_in .*= a0_surf
    ml_in[3] = (2.9 + i)*1e-10
    stheta, sphi, costheta, cosphi = sin(theta), sin(phi), cos(theta), cos(phi)
    
    ml_o = ml_in + [v*stheta*cosphi, v*stheta*sphi, v*costheta]
    ml_c = ml_in + [-w*stheta*cosphi, -w*stheta*sphi, -w*costheta]
    ml_bc = ml_in + [-bc*stheta*cosphi, -bc*stheta*sphi, -bc*costheta]
    
    out_attr = mol_surf_attr_stone(ml_o, ml_c, ml_bc, costheta)
    out_rep, out_disp = mol_surf_rep_stone(ml_o, ml_c,4)
    push!(en, [out_attr, out_disp, out_rep] .* joule2wn)
end

en = hcat(en...)
xval = 2.9 .+ collect(0:0.1:1.5)
using Plots
gr(leg=false)
plotly()
# Set default plot attributes
default(titlefontsize=20, linewidth=3, legendfontsize=10, 
        guidefont=font(16), tickfont=font(14), thickness_scaling=1.25, frame=:box, size = (1200, 800))
# plot(xval, (en[1,:] - en[2,:] + en[3,:]), xlabel = "z/Å", ylabel = "energy/cm-1", marker=:circle, 
#         ls=:solid, c=:black, fmt=:png, label = "total", legend=:bottomright, title = "θ=0°")
plot(xval, en[1,:], ls=:dot, c=:black, label = "attraction", xlabel = "z/Å", 
    ylabel = "energy/cm-1",legend=:bottomright, title = "θ=180°")
plot!(xval, en[3,:], ls=:dashdot, c=:black, label = "repulsion")
plot!(xval, -en[2,:], ls =:dash, c=:black, label = "dispersion")
#default()
plot!(xval, (en[1,:] - en[2,:] + en[3,:]), marker=:circle, 
ls=:solid, c=:black, label = "total",legend=:topright, inset = (1, bbox(0.02, 0.02, 0.5, 0.5, :top, :right)), subplot = 2)

all = []
for i in -20:10:200
    x[1] = i
    res2 = optimize(z_mono2, x, LBFGS(), Optim.Options(x_tol = 1e-8, f_tol = 1e-8, iterations = 50))
    p = [i, Optim.converged(res2), minimum(res2), 3.1+res2.minimizer[5]]
    println(p)
    push!(all, p)
end

all = hcat(all...)'

plot(collect(-20:10:200), all[:,3])


site_surface_interaction(com, ml_in, bl, theta1, phi1)

function z_mono(x)
    theta1 = x[1]*pi/180
    phi1 = x[2]*pi/180
    ml_in = zeros(Float64,3)
    com = zeros(Float64,3)
    ml_in[1] = x[3]
    ml_in[2] = x[4]
    com[3] = z0
   
    ml_in[1] = ml_in[1] - round(ml_in[1])
    ml_in[2] = ml_in[2] - round(ml_in[2])
    
    bl = v+w
    out = site_surface_interaction(com, ml_in, bl, theta1, phi1)
    return out[1]*joule2wn
end


pot = []
coord = []
for z in 2.4:0.1:4.2
    print("\n",z,"\t")
    global z0 = z*1e-10/a0_surf
    local x = zeros(Float64,4)
    x[1] = 32
    res=z_mono(x)
    print(res,"\t")
    push!(pot, res)
    
    #push!(pot, minimum(res))
    #ush!(coord, res.minimizer)
end


zval = collect(3:0.1:4.2)
cord = hcat(coord...)'
cord[:,3:4] = cord[:,3:4].*a0_surf/1e-10
out_val = hcat(zval, pot, cord)

# out_file = open("monolayer_z_pes-singlemol.txt", "w")
# using DelimitedFiles
# print(out_file,"z/Å","\t","min energy/cm-1\t","theta\t","phi\t x/Å\t y/Å \n")
# writedlm(out_file,out_val,"\t")
# close(out_file)


x = [180.0, 0, 0, 0, 1.3]
val = []
zval = collect(0.0:0.1:1.4)

val = [z_mono2([0.0, 0.0, 0.0, 0.0, i]) for i in zval]

trace = scatter(x = zval .+ 2.6, y=val,marker_color="black",mode="markers+lines")
plot(trace)

#site_surface_interaction([0,0,0.0], [0,0,3.5/3.99],v+w, 15*pi/180, 0.0) .*joule2wn

x = [0, 0, 0, 0, 0.35]
z_mono2(x)

site_surface_interaction([0.0,0.0,0.0], [0.0,0.0,0.79],   v+w, pi*1, 0.0) .*joule2wn
t = pi

ml_in = [0.0, 0.0, 3.3516e-10]
(mol_surf_rep(ml_in, t, 0.0, 0, 3))*joule2wn
(mol_surf_rep(ml_in, t, 0.0, 1, 3))*joule2wn
