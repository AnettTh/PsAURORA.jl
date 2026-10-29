using AURORA
using CairoMakie


## Make magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = DipoleMagneticField()

## Define the particles
μ = -cos(deg2rad(5))
L = 6.5
r0 = Cartesian(L*RE, 0.0, ϕ_eq)

E_grid = range(1e3, 40e3, length=100)
particles = [
    ParticleState(
        E,
        μ,
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=true,
    )
    for E in E_grid
]

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=500)
desnity = DentonDensity(magnetic_field, L, ϕ_eq)
plasma = PlasmaState(λ_grid, density)

## Define the wave
Ω_e_eq = minimum(plasma.Ω_e)
ωs = range(Ω_e_eq*0.1, Ω_e_eq*0.5, length=10)

modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

## Calculate TOF
TOF_matrix = zeros(length(E_grid), length(ωs))  # [n_E × n_ω]
TOF_simple = zeros(length(E_grid), length(ωs))


for (i, p) in enumerate(particles)
    TOF_matrix[i, :] = WPI_TOF(ωs, p, plasma; wave=wave, field_dependent=true )
    TOF_simple[i, :] = WPI_TOF(ωs, p, plasma; wave=wave, field_dependent=false)
end

## Make Figure
fig = Figure(size=(900,500))
ax = Axis(
    fig[1,1],
    xlabel="time-of-flight [s]",
    ylabel="E [keV]",
    title="Dipole approximation with pitch-angle α=$(round(rad2deg(acos(abs(μ)))))°",
)

#ax.yticks = [1, 10, 100, 1000]
#ax.ytickformat = values -> ["$(Int(v))" for v in values]

ω_norm = (ωs .- minimum(ωs)) ./ (maximum(ωs) - minimum(ωs))


#colors = [:blue, :red, :green, :orange, :purple]  # one per ω
colors = cgrad(:turbo, length(ωs), categorical=true)

for (i, ω) in enumerate(ωs)
    lines!(ax, TOF_matrix[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:solid,
        label="ω = $(round(ω * 2π/1e3, digits=1)) kHz"
    )
    lines!(ax, TOF_simple[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:dash
    )
end
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:solid,  label="Field-dependent")
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:dash,   label="Field-independent")

Legend(fig[1,2], ax)

##
#save("src/WPI/test_scripts/field_independent_model/figures/several_energies_$(round(rad2deg(acos(abs(μ))))).png", fig)
