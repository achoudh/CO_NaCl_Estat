include("initialization.jl")

σ = eigenvals ./ nmols_ml
esu(l, i) = eigenvecs[i, l] * (μ01 .* eu[i,:]) #each coloume is eigen vector l
pl = [sum([esu(l, i) for i in 1:nmols_ml]) for l in 1:nmols_ml]

μEpda = [(0.5*((pl[l][1])^2 + (pl[l][2])^2)*Tx + (pl[l][3])^2 *Tz) for l in 1:nmols_ml]
μEsda = [0.5*(pl[l][1]^2 + pl[l][2]^2)*Ty for l in 1:nmols_ml]

σpda = [unit2*σ[n1]*μEpda[n1] for n1 in 1:nmols_ml]
σsda = [unit2*σ[n1]*μEsda[n1] for n1 in 1:nmols_ml]

#using Plots
# plot(crossdataPda, crossdataSda, 
#      xlabel = "Frequency/cm⁻¹", 
#      ylabel = "Cross section/a.u.", 
#      label = ["Absorption" "Emission"], 
#      size = (600, 400), 
#      legend = :topleft, 
#      grid = :none)

     # Line-shape function (gaussian)
#plot(eigenvals, σpda)
# Domain Average

# Wavenumber array

range = 20  # "Wavenumbers"
dtponts = 200

step = 2 * (range / dtponts)

νk = collect(ν0[2] - range :step:ν0[1] + range)

# Function

# Line width

# const Δν = 1.1  # "Wavenumbers" # FWHM
# sσ(Δν) = Δν / (2log(4))

sσ(Δν) = Δν / (2sqrt(2log(2)))

gssn(x, ν, Δν) = exp(-(x - ν)^2 / (2 * sσ(Δν)^2)) / (sσ(Δν) * sqrt(2π))

ip(x) = sum([σpda[n] * gssn(x, eigenvals[n], 1.15) for n in 1:nmols_ml])

is(x) = sum([σsda[n] * gssn(x, eigenvals[n], 1.1) for n in 1:nmols_ml])

# Single domain

# ip1(x) = sum([σpda[n] * gssn(x, eigenvals[n], 1.1) for n in 1:nmols_ml])

# is1(x) = sum([σsda[n] * gssn(x, eigenvals[n], 1.1) for n in 1:nmols_ml])

# iip1 = hcat(νk, ip1.(νk))

# iis1 = hcat(νk, is1.(νk))

# # define iip1 and iis1 as arrays

# create the plot

#Create the traces for each spectrum
trace1 = scatter(x=νk, y=ip.(νk), mode="lines", name="p-pol")
trace2 = scatter(x=νk, y=is.(νk), mode="lines", name="s-pol")

#Create the plot layout
layout = Layout(
    title="Exciton result",
    xaxis_title="Wavelength",
    yaxis_title="Intensity",
    bgcolor="white"
)

# Combine the traces and layout into a plot object
plot_data = [trace1, trace2]
plot1 = plot(plot_data, layout)

# Display the plot in the notebook
display(plot1)
