function potco(r::Float64)
    # 2020/4/21  C-O AE-CCSD(T)(short range) MRCI(long range) / ACVTZ Energy in cm-1
    # updated: 2020/6/25, dv/dr
    ns::Int64 = 8
    w11::Vector{Float64} = [2.5731334302226245, -4.1635475296365723, 6.2652688635816842, -8.1228824615517947, 2.5634824563114837e1, -1.6643878666848999, 1.2863481593265147e1, -5.1395051186435685]
    b11::Vector{Float64} = [-1.6421783363943476, 1.5364774081764951, -1.1338455332512607, -1.4448051810696709, 7.5217573991947644, -1.4005229119446290, 1.1053854378210930e1, -5.9299626180269485]
    w21::Vector{Float64} = [-8.8664820626030757e-03, 7.8571245473773067e-03, -3.6411047563733342e-03, -4.0358215533209145e-03, 9.6640587889626451e-04, -1.4325782866595651, 1.2002568907875554e-02, 8.3983298757280007]
    b21::Float64 = 6.8970075393140338
    ra::Float64 = 1.4
    rb::Float64 = 7.0
    va::Float64 = 9.8087308941410573e-02
    vb::Float64 = 1.9558422718340193e+05

    v::Float64 = 0.0
    dv::Float64 = 0.0
    y1j::Float64 = 0.0

    x::Float64 = 2.0 * (r - ra) / (rb - ra) - 1.0
    for j in 1:ns
        y1j = tanh(b11[j] + w11[j] * x)
        v += w21[j] * y1j
        dv += w21[j] * (1.0 - y1j^2) * w11[j]
    end
    v += b21
    v = ((v + 1.0) * (vb - va) / 2.0) + va
    v += 0.560096850315234
    # dv *= (vb - va) / (rb - ra)
    return v #, dv
end

function nsimx(nt::Int64,r::Array{Float64,2},v::Array{Float64,1},idv::Int64,dv::Array{Float64,2},
    n0::Int64,n1::Int64,n2::Int64,rg::Array{Float64,2},w1::Array{Float64,2},b1::Array{Float64,1},
    w2::Array{Float64,2},b2::Array{Float64,1},w3::Array{Float64,1},b3::Float64,vg::Array{Float64,1})

    x0::Array{Float64,2} = r
    rgg = zeros(Float64, n0)
    for i in 1:n0
    rgg[i] = rg[2,i] - rg[1,i]
    x0[i,:] = (x0[i,:] .- rg[1,i]) ./ rgg[i] .* 2 .- 1
    end

    x1::Matrix{Float64} = fill(b1, (nt,n1))
    x1 .= x0' * w1 .+ x1

    x1 = dtanh.(x1)

    x2::Matrix{Float64} = fill(b2, (nt,n2))
    x2 .= x1 * w2 .+ x2

    x2 = dtanh.(x2)

    v::Vector{Float64} .= fill(b3,nt)
    v .= x2 * w3 .+ v

    v .= (v .+ 1) ./ 2 .* (vg[2] - vg[1]) .+ vg[1]

    if idv != 1
    return
    end

    dv::Matrix{Float64} .= 0.0
    for n in 1:nt
    for i in 1:n0
    for k in 1:n2
        tmp::Float64 = 0.0
        for j in 1:n1
            tmp += w2[j,k] * w1[i,j] * (1 - x1[n,j]^2)
        end
        dv[i,n] += w3[k] * tmp * (1 - x2[n,k]^2)
    end
    dv[i,n] *= (vg[2] / rgg[i])
    end
    end
end


function nsim(r0::Vector{Float64}, idv::Int, n0::Int, n1::Int, n2::Int, rg::Matrix{Float64}, w1::Matrix{Float64},
    b1::Vector{Float64}, w2::Matrix{Float64}, b2::Vector{Float64}, w3::Vector{Float64}, b3::Float64, vg::Vector{Float64})
    # ---- simulate a neural network with two hidden layers
    # ----      n0-n1-n2-1
    # ---- blas routines used in this subroutine: dgemv, ddot
    r::Vector{Float64} = copy(r0)
    v::Float64 = 0.0
    rgg = zeros(Float64, n0)
    vgg::Float64 = vg[2] - vg[1]
    
    # mapminmax [-1,1]
    for i in 1:n0
        rgg[i] = rg[2,i] - rg[1,i]
        r[i] = 2.0 * (r[i] - rg[1,i]) / rgg[i] - 1.0
    end
    
    # 1st layer
    rt1::Vector{Float64} = b1 .+ transpose(w1) * r
    ax::Vector{Float64} = tanh.(rt1)
    # 2nd layer
    rt2::Vector{Float64} = b2 .+ transpose(w2) * ax
    bx::Vector{Float64} = tanh.(rt2)

    # output layer
    v = b3 + dot(w3, bx)
    # reverse map
    v = vgg * (v[1] + 1.0) / 2.0 + vg[1]

    return v
end

# dtanh(x::Float64) = tanh(x) * (1.0 - tanh(x))
# dtanh(x::Vector{Float64}) = tanh.(x) .* (1.0 .- tanh.(x))

function nnfit_ococ(ndim, ntot, r0, idv)

    vx = zeros(Float64, ntot)
    # dv = zeros(Float64, ndim, ntot)

    if ndim != s0_nn
        error("ndim ≠ s0_nn")
    end

    vx = nsim(r0[:, 1], idv, s0_nn, s1_nn, s2_nn, rg_nn[:,:], 
              w1_nn[:,:], b1_nn[:], w2_nn[:,:], b2_nn[:], w3_nn[:], b3_nn, vg_nn)
   
    return vx
end


function mfi_ococ(r::Vector{Float64})# , idv::Int64)
    # molecular based fundamental invariants
    g::Vector{Float64} = r #exp.(-0.3 .* r)
    p = zeros(Float64, 7)
    # dpdr = zeros(7, 6)
    p[1] = g[1] + g[6]
    p[2] = g[3] + g[4]
    p[3] = g[1]^2 + g[6]^2
    p[4] = g[3]^2 + g[4]^2
    p[5] = g[1]*g[3] + g[4]*g[6]
    p[6] = g[2]
    p[7] = g[5]

    p[3] = sqrt(p[3])
    p[4] = sqrt(p[4])
    p[5] = sqrt(p[5])

    return p 
end

function pes_ococ(r::Vector{Float64})
    #pi = acos(-1.0)
    np::Int64 = 7
    p = zeros(Float64, np)
    # dv = zeros(Float64, np)
    vp = zeros(Float64, 1)

    z::Float64 = (r[2]+r[3]+r[4]+r[5])/4.0

    p = mfi_ococ(r)
    
    vp = nnfit_ococ(np, 1, p, 0)

    v::Float64 = vp[1]

    if z > 20.0
        v = 0.0
    elseif z > 18.0 # 18-20
        v = v * (0.5 - 0.5*sin((z-19.0)*pi*0.5))
    end

    # v = vp[1]
    return v #, dv
end

function co_co_NNpes(R12::Vector{Float64}, phi1::Float64, theta1::Float64, phi2::Float64, theta2::Float64)
    # Determines the electrostatic interaction energy between two CO molecules

    # local axis unit vectors
    stheta1::Float64, sphi1::Float64, costheta1::Float64, cosphi1::Float64 = sin(theta1), sin(phi1), cos(theta1), cos(phi1)
    stheta2::Float64, sphi2::Float64, costheta2::Float64, cosphi2::Float64 = sin(theta2), sin(phi2), cos(theta2), cos(phi2)
    
    e1z::Vector{Float64} = [stheta1 * cosphi1, stheta1 * sphi1, costheta1]
    e2z::Vector{Float64} = [stheta2 * cosphi2, stheta2 * sphi2, costheta2]

    # vv = 1.128e-10 * 12/28
    # ww = 1.128e-10 * 16/28
    
    o1p::Vector{Float64} = (v .* e1z) #[ v*sin(theta1)*cos(phi1),  v*sin(theta1)*sin(phi1),  v*cos(theta1)]
    o2p::Vector{Float64} = (v .* e2z) #[ v*sin(theta2)*cos(phi2),  v*sin(theta2)*sin(phi2),  v*cos(theta2)]
    c1p::Vector{Float64} = (-w .* e1z) #[-w*sin(theta1)*cos(phi1), -w*sin(theta1)*sin(phi1), -w*cos(theta1)]
    c2p::Vector{Float64} = (-w .* e2z) #[-w*sin(theta2)*cos(phi2), -w*sin(theta2)*sin(phi2), -w*cos(theta2)]
    
    # Distance between the atoms from two molecules
    o1c1::Float64 = norm(o1p-c1p)
    o2c2::Float64 = norm(o2p-c2p)
    o1o2::Float64 = norm(R12 .+ (o1p.-o2p))
    o1c2::Float64 = norm(R12 .+ (o1p.-c2p))
    c1c2::Float64 = norm(R12 .+ (c1p.-c2p))
    o2c1::Float64 = norm(R12 .+ (c1p.-o2p))
    
    r::Vector{Float64} = [o1c1, o1o2, o1c2, o2c1, c1c2, o2c2] / bohr_radius 

    V_CO_CO_nu::Float64 = pes_ococ(r) # + potco(r[1])[1] + potco(r[6])[1]

    return V_CO_CO_nu
end


function co_co_int_nu_guo(r1, r2, r0, theta1, theta2, phi)
    #    V_CO_CO, Vqq, Vqmu, Vmumu, VmuQu, VmuO, VQuQu, VQuO, VOO, V_rep)
    # Determines the electrostatic interaction energy between two CO molecules
    θ1 = theta1 * pi/180 
    θ2 = theta2 * pi/180 
    ϕ  = phi * pi/180

    rcm1 = [0.0, 0.0, 0.0]
    rcm2 = [0.0, 0.0, r0]
    
    mc, mo = 12.011, 15.999
    M = mc+mo
    
    d_o = (mc / M)
    d_c = (mo / M)

    o1p = [0.0,  d_o * r1 * cos(θ1-pi/2.0), - d_o * r1 * sin(θ1-pi/2.0)]
    c1p = [0.0, - d_c * r2 * cos(θ1-pi/2.0),  d_c * r2 * sin(θ1-pi/2.0)]
    println(o1p)
    println(c1p)

    o2p = [   d_o * r1 * sin(θ2) * sin(ϕ),   d_o * r1 * sin(θ2) * cos(ϕ), r0 + d_o * r1 * cos(θ2)]
    c2p = [ - d_c * r2 * sin(θ2) * sin(ϕ), - d_c * r2 * sin(θ2) * cos(ϕ), r0 - d_c * r2 * cos(θ2)]
    println(o2p)
    println(c2p)
    
    println(acos( dot(rcm2, (o1p-c1p)) / (r0*r1))  *180 /pi)
    println(acos( dot(rcm2, (o2p-c2p)) / (r0*r2))  *180 /pi)

    # Distance between the atoms from two molecules
    o1c1 = norm(o1p-c1p)
    o2c2 = norm(o2p-c2p)
    o1o2 = norm(o1p-o2p)
    o1c2 = norm(o1p-c2p)
    c1c2 = norm(c1p-c2p)
    o2c1 = norm(c1p-o2p)

    r = [o1c1, o1o2, o1c2, o2c1, c1c2, o2c2] *1e-10 / bohr_radius

    println(r)

    # pic = scatter3d([o1p[1], o2p[1]], [o1p[2], o2p[2]], [o1p[3], o2p[3]], c=:red, ms=8, xlabel="x", ylabel="y", zlabel="z")
    # pic = scatter3d!([c1p[1], c2p[1]], [c1p[2], c2p[2]], [c1p[3], c2p[3]], c=:black, ms=9, ticks = nothing)
    # pic = scatter3d!([rcm1[1], rcm2[1]], [rcm1[2], rcm2[2]], [rcm1[3], rcm2[3]], c=:blue, ms=3, ticks = nothing)
    # #scatter3d!(ml_c[:,1], ml_c[:,2], ml_c[:,3], c=:black, ms=9, ticks = nothing)

    # xlims!(-1,1)
    # #zlims!(-4,4)
    # ylims!(-1,1)

    # display(pic)

    V_CO_CO_nu = pes_ococ(r) #+ potco(r[1]) + potco(r[6])

    return V_CO_CO_nu
end


# @time co_co_int_nu_guo(1.1282, 1.1282, 4.3335, 134.95, 45.05, 180.0)
# co_co_int_nu_guo(1.1282, 1.1282, 3.6181, 65.20, 114.80, 180.0)

# r = [2.131989014741522, 9.568805089619767, 7.977567964771465, 7.977567964771468, 6.694169458560173, 2.131989014741522]

# pes_ococ(r)
# pip_ococ(r)