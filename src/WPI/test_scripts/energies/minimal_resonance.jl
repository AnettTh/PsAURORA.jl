using AURORA
using CairoMakie

# NOTE: This does not look right
## Make magnetic field
ϕ = deg2rad(120.0)
tsyg_field = TsyganenkoMagneticField(ϕ)
dipole_field = DipoleMagneticField()


## Make a particle
μ = -cos(deg2rad(3))
L = 6.0
r0 = Cartesian(L*RE, 0.0, 0.0)
E = 30e3

particle_tsyg = ParticleState(E, μ, r0, tsyg_field; relativistic=true)
particle_dipole = ParticleState(E, μ, r0, dipole_field; relativistic=true)

## Make a plasma
λ_grid = range(0.0, deg2rad(50), length=500)
plasma_tsyg = PlasmaState(λ_grid, ϕ, DentonDensity(tsyg_field, L, ϕ), tsyg_field, Float64(L))
plasma_dipole = PlasmaState(λ_grid, ϕ, DentonDensity(dipole_field, L, ϕ), dipole_field, Float64(L))


## Make a wave
ω_tsyg = range(plasma_tsyg.Ω_e[1]*0.1, plasma_tsyg.Ω_e[1]*0.4, length=100)
ω_dipole = range(plasma_dipole.Ω_e[1]*0.1, plasma_dipole.Ω_e[1]*0.4, length=100)


## Find minimal resonance energy
E_min_tsyg = minimal_resonance_energy(ω_tsyg, plasma_tsyg, particle_tsyg)
E_min_dipole = minimal_resonance_energy(ω_dipole, plasma_dipole, particle_dipole)

## Make the figure
fig = Figure()
ax = Axis(fig[1, 1], xlabel="ω [rad/s]", ylabel="Eₘᵢₙ [keV]", title="Threshold for cyclotron resonance, ϕ = $(round(rad2deg(ϕ), digits=1))")

lines!(ax, ω_tsyg, E_min_tsyg / eV_in_J / 1e3, label="Tsyganenko")
lines!(ax, ω_dipole, E_min_dipole / eV_in_J / 1e3, label="Dipole")

axislegend(ax; position=:rt)


#save("src/WPI/test_scripts/energies/Emin_phi_$(round(rad2deg(ϕ), digits=1)).png", fig)
