function monolayer(nx, ny, thetamono, phimono)
    n = 1
    com::Matrix{Float64} = zeros(nmols_uc*nx*ny, 3)
    theta::Vector{Float64} = zeros(nmols_uc*nx*ny)
    phi::Vector{Float64} = zeros(nmols_uc*nx*ny)

    for i in 0:nx-1
        for j in 0:ny-1
            com[n:n+nmols_uc-1,1] .= unitcell[:,1] .+ 2*i
            com[n:n+nmols_uc-1,2] .= unitcell[:,2] .+ 2*j
            com[n:n+nmols_uc-1,3] .= unitcell[:,3]
            n += nmols_uc
        end
    end

    for i in 1:nmols_uc
        theta[i] = thetamono
        phi[i] = phimono[i]
    end

    i = nx*ny
    j = nmols_uc*i
    theta = vec(repeat(theta[1:nmols_uc], outer=(1,i)))
    phi = vec(repeat(phi[1:nmols_uc], outer=(1,i)))
    #phi = reshape(repeat(phi[1:nmols_uc], outer=(1,i)), j)

    return com, theta, phi
end

#, CSV, DataFrames

function display_structure(comml, θ_ml, ϕ_ml)
    
    cs = 20
    cls = 20
    nac = 10

    ll = nx*2 - 1
    mm = ny*2 - 1

    posna= [[i, j, 0] for i in 0:ll for j in 0:mm]

    poscl = [[i+0.5, j+0.5, 0] for i in 0:ll-1 for j in 0:mm-1]

    #na_x = [posna[i][1]  for i in eachindex(posna)]
    na_r = hcat(posna...)
    cl_r = hcat(poscl...)

    # na = scatter3d(x=na_r[1,:], y=na_r[2,:], z=na_r[3,:], mode="markers", marker_size=nac, marker_color ="blue")
    # cl = scatter3d(x=cl_r[1,:], y=cl_r[2,:], z=cl_r[3,:], mode="markers", marker_size=cls, color =:green)

    vv = 0.483/a0_surf
    ww = 0.645/a0_surf

    mlcxyz = zeros(Float64,nmols_ml,3)
    for i in 1:nmols_ml
        mlcxyz[i, 1] = comml[i, 1] - ww * sin(θ_ml[i]) * cos(ϕ_ml[i])
        mlcxyz[i, 2] = comml[i, 2] - ww * sin(θ_ml[i]) * sin(ϕ_ml[i])
        mlcxyz[i, 3] = comml[i, 3] - ww * cos(θ_ml[i]) + 1
    end
    
    mloxyz = zeros(Float64,nmols_ml,3)
    for i in 1:nmols_ml
        mloxyz[i, 1] = comml[i, 1] + vv * sin(θ_ml[i]) * cos(ϕ_ml[i])
        mloxyz[i, 2] = comml[i, 2] + vv * sin(θ_ml[i]) * sin(ϕ_ml[i])
        mloxyz[i, 3] = comml[i, 3] + vv * cos(θ_ml[i]) +1
    end

    c_ml = scatter3d(x=mlcxyz[:,1], y=mlcxyz[:,2], z=mlcxyz[:,3], mode="markers", marker_size=cs, marker_color="black")
    o_ml = scatter3d(x=mloxyz[:,1], y=mloxyz[:,2], z=mloxyz[:,3], mode="markers", marker_size=cs, marker_color="red")

    layout = Layout(scene=attr(xaxis_title="X", yaxis_title="Y", zaxis_title="Z", zaxis_range=(0,7) ),title="Structure")#zaxis=attr(range=(0,100)),

    p2 = display(plot([na, cl, c_ml, o_ml], layout))
    return p2
end

# define fij
function fij(i, j)
    rn = com[i,:] - com[j,:]
    rn[1] = rn[1] - grid * round(rn[1] / grid)
    rn[2] = rn[2] - grid * round(rn[2] / grid)
    rn[3] = 0
    rn = rn .* a0_surf
    r = norm(rn)
    #print(rn)
    en = rn ./ r
    return 1 / r^3 * (dot(eu[i,:], eu[j,:]) - 3 * dot(en, eu[i,:]) * dot(en, eu[j,:]))
end

# define sum1
function sum1(n1)
    s = 0
    for j = 1:nmols_ml
        if j != n1
            #print("$j\t")
            s += fij(j, n1)
            #println(fij(j,n1))
        end
    end
    return s
end

