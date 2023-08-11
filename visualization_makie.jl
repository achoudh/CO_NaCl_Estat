# using Plots
# #gr()
# plotlyjs()

# # Set default plot attributes
# default(titlefontsize=12, labelfontsize=10, linewidth=2, legendfontsize=8, 
#         guidefont=font(10), tickfont=font(10), thickness_scaling=1.15, size = (700, 400)) #, frame=:box

# scene = meshscatter(x, y, z, markersize = 10, color = :white)
# scene
# display(scene)

size_reduce = 0.02

function nacl_show()
        ll = nx*2 - 1
        mm = ny*2 - 1

        posna = hcat([[i, j, 0] for i in 0:ll for j in 0:mm]...)
        poscl = hcat([[i+0.5, j+0.5, 0] for i in 0:ll-1 for j in 0:mm-1]...)

        return posna, poscl
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
        
return ml_c, ml_o, ol_c, ol_o
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

function write_to_file(file_path, data)
        # Step 1: Open the file in write mode
        file = open(file_path, "w")
        println(file, "number of iterations")
        println(file, size(data[1],1)-1)

        println(file,"\nBest Energies from each iterations")
        writedlm(file, [data[1]])
        
        println(file,"\nBest states from each iterations, state 1 in initial")
        writedlm(file, size(data[2][1],1))
        writedlm(file, data[2])

        println(file,"\nEnergies at accepted steps")
        writedlm(file, size(data[3],1))
        writedlm(file, [data[3]])

        close(file)
        
end

function show_figure(x::Vector{Float64}, lattice_ml::Matrix{Float64}, lattice_ol::Matrix{Float64}, figtitle::String, plot_ol::Int64 = 1)
        fig0 = Figure(resolution=(600,800))

        ax_sp = Axis(fig0[1:2, 1],
                title  = figtitle,
                titlesize = 22,
                titlefont = "Times New Roman",
                subtitle = L"IR-Spectra (domain averaged)$ $",
                xlabel = L"Frequnecy/cm$^{-1}$",
                ylabel = L"Intentsity/a.u.$ $",
                xlabelsize = 18, 
                ylabelsize = 18, 
                xticklabelfont = "Times New Roman",
                yticklabelfont = "Times New Roman")

        ax_st = LScene(fig0[3:5, 1], show_axis=false) #, title = L"Structure$ $") ,limits = Rect(0,0,0,1,1,1)

        ipda, isda, ip, is = ir_spectra(νk, x, lattice_ml, Δν)
        lines!(ax_sp, νk, ipda, color=:black, label=L"p-pol$ $")
        lines!(ax_sp, νk, isda, color=:blue,  label=L"s-pol$ $")
        axislegend(ax_sp, titlefont = "Times New Roman")

        ml_C, ml_O, ol_C, ol_O  = structure_unitmono(x, lattice_ml, lattice_ol)
        meshscatter!(ax_st, ml_C[:,1], ml_C[:,2], ml_C[:,3], markersize = 0.2, color=:black)#, limits = Rect(-1, -1, -1, 5, 5, 1))
        meshscatter!(ax_st, ml_O[:,1], ml_O[:,2], ml_O[:,3], markersize = 0.2, color=:red)
        if plot_ol==1
                meshscatter!(ax_st, ol_C[:,1], ol_C[:,2], ol_C[:,3], markersize = 0.2, color=(:black, 0.3))#, transparency = true)
                meshscatter!(ax_st, ol_O[:,1], ol_O[:,2], ol_O[:,3], markersize = 0.2, color=(:red, 0.3))#, transparency = true)
        end

        r_Na, r_Cl = nacl_show()
        meshscatter!(ax_st, r_Na[1,:], r_Na[2,:], r_Na[3,:], markersize = 0.1, color=:gray)
        meshscatter!(ax_st, r_Cl[1,:], r_Cl[2,:], r_Cl[3,:], markersize = 0.3, color=:green)
        return fig0
end

function set_axis(fig, xlabel, ylabel, ratio)
        ax = Axis(fig,
                titlesize = 22*ratio,
                titlefont = "Times New Roman",
                xlabel = xlabel,
                ylabel = ylabel,
                xlabelsize = 18*ratio,
                xlabelfont = "Times New Roman",
                ylabelfont = "Times New Roman",
                ylabelsize = 18*ratio, 
                xticklabelfont = "Times New Roman",
                yticklabelfont = "Times New Roman",
                xticklabelsize = 16*ratio,
                yticklabelsize = 16*ratio)

                return ax
end