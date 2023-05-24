using Plots
gr()
# Set default plot attributes
default(titlefontsize=20, labelfontsize=16, linewidth=2, legendfontsize=14, 
        guidefont=font(16), tickfont=font(14), thickness_scaling=1.15, frame=:box, size = (900, 600))

function nacl_show(r_Cl,r_Na)

        plotly()

        ll = nx*2 - 1
        mm = ny*2 - 1

        posna = hcat([[i, j, 0] for i in 0:ll for j in 0:mm]...)
        poscl = hcat([[i+0.5, j+0.5, 0] for i in 0:ll-1 for j in 0:mm-1]...)

        surface = scatter3d(posna[1,:], posna[2,:], posna[3,:], ms = r_Na, c=:violet, label = "Na")
                 scatter3d!(poscl[1,:], poscl[2,:], poscl[3,:], ms = r_Cl, c=:green,  label = "Cl", alpha=0.8)
        return surface
end


function display_structure_unitmono(x, lattice_ml)
        
        vv, ww = [v, w] ./a0_surf
        
        ml_o = zeros(Float64,nmols_ml,3)
        ml_c = zeros(Float64,nmols_ml,3)
        
        for i in 1:nmols_ml

                stheta, sphi, costheta, cosphi = sin(x[i]), sin(x[i+nmols_ml]), cos(x[i]), cos(x[i+nmols_ml])
                rvec = [x[i+2*nmols_ml], x[i+3*nmols_ml], x[i+4*nmols_ml]] + lattice_ml[i,:]
            
                ml_o[i,:]  = rvec + [vv*stheta*cosphi, vv*stheta*sphi, vv*costheta]
                ml_c[i,:]  = rvec + [-ww*stheta*cosphi, -ww*stheta*sphi, -ww*costheta]

        end
        
        co_nacl = nacl_show(20, 10)
#        display(nacl_show(20, 10))
        scatter3d!(co_nacl,ml_c[:,1],ml_c[:,2],ml_c[:,3],ms=9,c=:black,camera=(0,90,),label="C")
        scatter3d!(co_nacl,ml_o[:,1],ml_o[:,2],ml_o[:,3],ms=8,c=:red, ticks=nothing, label = "O")
        zlims!(co_nacl,0,1)
return display(co_nacl)
end
        