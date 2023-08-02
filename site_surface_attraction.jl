include("co_co_interactions.jl") #This is to get the multipole_R_dependence


function surf_pot_arnab(com::Vector{Float64}, costheta::Float64, mom::Vector{Float64})::NTuple{4, Float64}

    x, y, z = com/a0_NaCl
    ret = zeros(Float64, 5)
    for i::Int64 in 1:nlm
        l::Int64, m::Int64 = lodd[i], modd[i]
        temp = ret_c[i] * exp(-2*pi*z*lm[i]) * cos(2*pi*((l*x)+(m*y)-(l+m)/4))
        ret[1] += temp
        for j in 2:4  # upto 4th moment
           temp *= lma[i]
           ret[j] -= temp 
        end
    end
    
    return mom[1]*ret[1], -costheta*mom[2]*ret[2], -(3*(costheta)^2-1)*mom[3]*ret[3]/4,
         -(5*(costheta)^3-3*costheta)*mom[4]*ret[4]/12
end
    

function mol_surf_attr_arnab(com::Vector{Float64}, costheta::Float64)::Float64
    op::Vector{Float64} = rot*com
    mom = collect(multipole_R_dependence(v+w))
    return sum(surf_pot_arnab(op, costheta, mom))
end
