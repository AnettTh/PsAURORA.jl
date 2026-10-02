using AURORA
using CairoMakie

# NOTE: This does not look right
## Make magnetic field
ϕ = 0.0
magnetic_field = TsyganenkoMagneticField(ϕ)

## Make a particle
μ = -cos(deg2rad(3))
L = 6.0
r0 = Cartesian(L*RE, 0.0, 0.0)
E = 30e3

particle = ParticleState(E, μ, r0, magnetic_field; relativistic=true)

## Make a plasma
λ_grid = range(0.0, deg2rad(50), length=500)
plasma = PlasmaState(λ_grid, ϕ, DentonDensity(magnetic_field, L, ϕ), magnetic_field, Float64(L))

## Make a wave
ω = range(plasma.Ω_e[1]*0.1, plasma.Ω_e[1]*0.4, length=100)


## Find minimal resonance energy
E_min = minimal_resonance_energy(ω, plasma, particle)   # NOTE: This might be in J? Check equation in paper.

## Make the figure
fig = Figure()
ax = Axis(fig[1, 1], xlabel="ω [rad/s]", ylabel="Eₘᵢₙ [eV]")

lines!(ax, ω, E_min)

fig
