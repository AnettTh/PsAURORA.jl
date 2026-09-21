using AURORA
using CairoMakie

# NOTE: For several energies
## Define the particles
μ = -cos(deg2rad(0))       # Almost field-aligned
L = 6.0
r0 = [L*RE, 0.0, 0.0]

E_grid = range(1e3, 10000e3, length=100)
particles = [ParticleState(E, μ, r0, dipole_field) for E in E_grid]

# Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
n_e0 = 1.8e7        # 18/cc, from Hsieh 2022

plasma = PlasmaParameters(λ_grid, n_e0, dipole_field, L)


# Define the wave
ωs = range(plasma.Ω_e[1]*0.2, plasma.Ω_e[1]*0.4, length=5)


## Calculate TOF
TOF_matrix = zeros(length(E_grid), length(ωs))  # [n_E × n_ω]
TOF_simple = zeros(length(E_grid), length(ωs))


# TODO: Figure out why only >100 keV precipitates??
for (i, p) in enumerate(particles)
    TOF_matrix[i, :] = WPI_TOF(ωs, p, plasma; field_dependent=true)
    TOF_simple[i, :]  = WPI_TOF(ωs, p, plasma; field_dependent=false)
end


## Make Figure
fig = Figure()
ax = Axis(fig[1,1], xlabel="time-of-flight [s]", ylabel="E [keV]", yscale=log10, title="Pitch-angle $(round(rad2deg(acos(abs(μ)))))")

ax.yticks = [1, 10, 100, 1000]
ax.ytickformat = values -> ["$(Int(v))" for v in values]

ω_norm = (ωs .- minimum(ωs)) ./ (maximum(ωs) - minimum(ωs))


colors = [:blue, :red, :green, :orange, :purple]  # one per ω

for (i, ω) in enumerate(ωs)
    lines!(ax, TOF_matrix[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:solid,
        label="ω = $(round(ω/1e3, digits=1)) kHz"
    )
    lines!(ax, TOF_simple[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:dash
    )
end
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:solid,  label="Field-dependent")
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:dash,   label="Field-independent")

axislegend(ax, position=:rb)

##
xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

save("src/WPI/test_scripts/several_energies_$(round(rad2deg(acos(abs(μ))))).png", fig)
