using AURORA
using CairoMakie


## Make magnetic field
ϕ_eq = deg2rad(240.0)
tsyg_field = TsyganenkoMagneticField(ϕ_eq)
dipole_field = DipoleMagneticField()


## Make a particle
μ = -cos(deg2rad(3))
L = 7.1
r0 = Cartesian(L*RE, 0.0, ϕ_eq)

# NOTE: This only uses the position of the particle, so the energy is irrelevant
particle_tsyg = ParticleState(1e4, μ, r0, tsyg_field, ϕ_eq; relativistic=true)
particle_dipole = ParticleState(1e4, μ, r0, dipole_field, ϕ_eq; relativistic=true)

## Make a plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=500)

density_tsyg = DentonDensity(tsyg_field, L, ϕ_eq)
density_dipole = DentonDensity(dipole_field, L, ϕ_eq)

plasma_tsyg = PlasmaState(λ_grid, density_tsyg)
plasma_dipole = PlasmaState(λ_grid, density_dipole)


## Make a wave
ω_tsyg = range(minimum(plasma_tsyg.Ω_e)*0.1, minimum(plasma_tsyg.Ω_e)*0.5, length=100)
ω_dipole = range(minimum(plasma_dipole.Ω_e)*0.1, minimum(plasma_dipole.Ω_e)*0.5, length=100)


## Find minimal resonance energy
E_min_tsyg = minimal_resonance_energy(ω_tsyg, plasma_tsyg, particle_tsyg)
E_min_dipole = minimal_resonance_energy(ω_dipole, plasma_dipole, particle_dipole)

## Make the figure
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="ω [rad/s]",
    ylabel="Eₘᵢₙ [keV]",
    title="Threshold for cyclotron resonance, ϕ = $(round(rad2deg(ϕ_eq), digits=1))"
)

lines!(ax, ω_tsyg, E_min_tsyg / eV_in_J / 1e3, label="Tsyganenko")
lines!(ax, ω_dipole, E_min_dipole / eV_in_J / 1e3, label="Dipole")

axislegend(ax; position=:rt)


#save("src/WPI/test_scripts/energies/Emin_phi_$(round(rad2deg(ϕ), digits=1)).png", fig)
