using AURORA
using AURORA; c₀
using CairoMakie

## Define magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
magnetic_field = DipoleMagneticField()

## Define the particle
μ = -cos(deg2rad(2))       # Almost field-aligned
L = 6.0
r0 = Cartesian(L*RE, 0.0, ϕ_eq)

E = 10e3
particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true)

## Define the plasma
λ_grid = range(deg2rad(-50), deg2rad(50), length=500)
lb_low, lb_high = 0.1, 0.5
density = DentonDensity(magnetic_field, L, ϕ_eq)
plasma = PlasmaState(
    λ_grid,
    density;
    lb_low=lb_low,
    lb_high=lb_high
    )


## Get group velocity
ω_grid=range(plasma.ω_lb[1], plasma.ω_lb[end], length=500)

# NOTE: I do not think this is correct anymore
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

## Calculate TOF
TOF_fd     = WPI_TOF(ω_grid, particle, plasma; wave, field_dependent=true)
TOF_simple = WPI_TOF(ω_grid, particle, plasma; wave, field_dependent=false)


## Make Figure
fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel="ω / Ωₑ",
    ylabel="TOF [s]",
    title="Dipole approximation with pitch-angle α=$(round(rad2deg(acos(abs(μ)))))°"
)

lines!(ax, [NaN], [NaN]; color=:transparent, label="E = $(round(E/1e3, digits=1)) keV")
lines!(ax, ω_grid / minimum(plasma.Ω_e), TOF_fd;     linestyle=:dash, label="Field-dependent")
lines!(ax, ω_grid / minimum(plasma.Ω_e), TOF_simple; linestyle=:solid,  label="Field-independent")

axislegend(ax, position=:lt)

##
save("src/WPI/test_scripts/field_independent_model/figures/several_ω_E$(round(E/1e3, digits=1)).png", fig)
