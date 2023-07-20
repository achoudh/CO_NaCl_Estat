using Plots
gr()
plotly()

# Set default plot attributes
default(titlefontsize=20, labelfontsize=16, linewidth=2, legendfontsize=14, 
        guidefont=font(16), tickfont=font(14), thickness_scaling=1.15, frame=:box, size = (900, 600))

function nacl_show(r_Cl,r_Na)

        ll = nx*2 - 1
        mm = ny*2 - 1

        posna = hcat([[i, j, 0] for i in 0:ll for j in 0:mm]...)
        poscl = hcat([[i+0.5, j+0.5, 0] for i in 0:ll-1 for j in 0:mm-1]...)

        surface = scatter3d(posna[1,:], posna[2,:], posna[3,:], ms = r_Na, c=:gray, label = "Na")
                 scatter3d!(poscl[1,:], poscl[2,:], poscl[3,:], ms = r_Cl, c=:green,  label = "Cl", alpha=0.8)
        return surface
end


function structure_unitmono(x, lattice_ml, lattice_ol)
        
        theta_ml = x[1+0*nmols_ml:1*nmols_ml]
        phi_ml =   x[1+1*nmols_ml:2*nmols_ml]
        δr_ml =    x[1+2*nmols_ml:5*nmols_ml]
        δz_ol =    x[ndofs_ml]
        theta_ol = x[1+ndofs_ml+0*nmols_ol2:ndofs_ml+1*nmols_ol2]
        phi_ol   = x[1+ndofs_ml+1*nmols_ol2:ndofs_ml+2*nmols_ol2]
        xy_ol    = x[1+ndofs_ml+2*nmols_ol2:ndofs_ml+4*nmols_ol2]
        
        vv, ww = [v, w] ./a0_surf
        
        ml_o = zeros(Float64,nmols_ml,3)
        ml_c = zeros(Float64,nmols_ml,3)
        ol_o = zeros(Float64,nmols_ol2,3)
        ol_c = zeros(Float64,nmols_ol2,3)
        
        for i in 1:nmols_ml
                stheta, sphi, costheta, cosphi = sin(theta_ml[i]), sin(phi_ml[i]), cos(theta_ml[i]), cos(phi_ml[i])
                rvec = [δr_ml[i], δr_ml[i+nmols_ml], δr_ml[i+2*nmols_ml] ] + lattice_ml[i,:]
                    
                ml_o[i,:]  = rvec + [vv*stheta*cosphi, vv*stheta*sphi, vv*costheta]
                ml_c[i,:]  = rvec + [-ww*stheta*cosphi, -ww*stheta*sphi, -ww*costheta]
        end

        for i in 1:nmols_ol2
            
                stheta, sphi, costheta, cosphi = sin(theta_ol[i]), sin(phi_ol[i]), cos(theta_ol[i]), cos(phi_ol[i])
                rvec = lattice_ol[i,:] + [xy_ol[i], xy_ol[i+nmols_ol2], δz_ol] 
            
                ol_o[i,:]  = rvec + [vv*stheta*cosphi, vv*stheta*sphi, vv*costheta]
                ol_c[i,:]  = rvec + [-ww*stheta*cosphi, -ww*stheta*sphi, -ww*costheta]

        end
        
        co_nacl = nacl_show(20, 10)
        scatter3d!(co_nacl,ml_c[:,1],ml_c[:,2],ml_c[:,3],ms=9,c=:black,camera=(0,90,),label="C")
        scatter3d!(co_nacl,ml_o[:,1],ml_o[:,2],ml_o[:,3],ms=8,c=:red, ticks=nothing, label = "O")
        scatter3d!(co_nacl,ol_c[:,1],ol_c[:,2],ol_c[:,3],ms=9,c=:black,camera=(0,90,),label="C",alpha = 0.3)
        scatter3d!(co_nacl,ol_o[:,1],ol_o[:,2],ol_o[:,3],ms=8,c=:red, ticks=nothing, label = "O",alpha = 0.3)
        # zlims!(co_nacl,0,10)

return co_nacl
end



using Printf
function show_params(x)
        
    θ_ml = x[1+0*nmols_ml:1*nmols_ml] /degrees #fill(38.0,4) * pi / 180.0
    ϕ_ml = x[1+1*nmols_ml:2*nmols_ml] /degrees
    ml_in = zeros(Float64, nmols_ml, 3)
    ml_in[:,1] = x[1+2*nmols_ml:3*nmols_ml]
    ml_in[:,2] = x[1+3*nmols_ml:4*nmols_ml]
    ml_in[:,3] = x[1+4*nmols_ml:5*nmols_ml] #fill(x[17],4)
    δz_ol = x[1+5*nmols_ml]

    out = hcat(θ_ml, ϕ_ml, ml_in[:,1].*a0_surf/1e-10, ml_in[:,2].*a0_surf/1e-10, (z_ml.+ml_in[:,3]).*a0_surf/1e-10)
    
    @printf("%-10s %-10s %-10s %-10s %-10s\n", "θ/°", "ϕ/°", "δx/Å", "δy/Å", "z/Å")
    for i in 1:nmols_ml
    @printf("%f %f  %f  %f  %f\n",out[i,1],out[i,2],out[i,3],out[i,4],out[i,5])
    end
    # return out
end