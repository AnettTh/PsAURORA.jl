using AURORA
using JLD2
using Dates
using CairoMakie

##
@load "src/WPI/data/11-25-2026/energy-time-position.jld2"

# Make E with the same shape as dx
E = collect(E_grid)
E_all = repeat(E, 1, size(dx_km, 2), size(dx_km, 3), size(dx_km, 4))
dx_all = fp_enu[1, :, :, :, :]

# Flatten and keep precipitated cases only
precipitated = findall(isfinite, vec(east))
dx = vec(dx_all)[precipitated]
E = vec(E_all)[precipitated]

# Filter and shift from time-of-flight to time-of-arrival-after-first-arrival
tof = vec(TOF)[precipitated]
t = tof .- minimum(tof)

##
fig = Figure(size=(800, 400))

Label(fig[0, 1:2], text="Energy as a function of east drift and time-of-arrival", font=:bold, fontsize=16)

ax1 = Axis(fig[1,1], xlabel="dx [km]", ylabel="E [keV]")
sc1 = scatter!(ax1, dx ./ 1e3, E ./ 1e3, markersize=4, alpha=0.5, color="#E69F00")

ax2 = Axis(fig[1,2], xlabel="t [s]")
sc2 = scatter!(ax2, t, E ./ 1e3, markersize=4, alpha=0.5, color="#E69F00")

linkyaxes!(ax1, ax2)
hideydecorations!(ax2, grid = false)

ylims!(ax1, 7, 40)
xlims!(ax1, 0, maximum(dx ./ 1e3))
xlims!(ax2, 0, maximum(t))

#save("src/WPI/test_scripts/input_flux/figures/11262026_DriftAndTimeVsEnergy.png", fig)
