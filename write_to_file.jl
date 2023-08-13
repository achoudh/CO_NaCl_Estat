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