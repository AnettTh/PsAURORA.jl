using AURORA
using AURORA; c₀
using CairoMakie

## Define the particle
μ = -cos(deg2rad(0))       # Almost field-aligned
L = 6.0
r0 = [L*RE, 0.0, 0.0]

E = 3e4
particle = ParticleState(E, μ, r0, dipole_field; relativistic=true)

## Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
n_e0 = 1.8e7        # 18/cc, from Hsieh 2022

plasma = PlasmaState(λ_grid, n_e0, dipole_field, L)


## Get group velocity
ω_grid=range(1e3, 4.74e3, length=500)

## Calculate TOF
TOF_fd     = WPI_TOF(ω_grid, particle, plasma; wave_launch_time=wave_chirp, field_dependent=true)
TOF_simple = WPI_TOF(ω_grid, particle, plasma; wave_launch_time=wave_chirp, field_dependent=false)


## Make Figure
fig = Figure()
ax = Axis(fig[1,1], xlabel="ω [kHz]", ylabel="time-of-flight [s]", title="Pitch-angle $(round(rad2deg(acos(abs(μ)))))")

lines!(ax, (ω_grid.*2π) ./ 1e3, TOF_fd;     linestyle=:dash, label="Field-dependent")
lines!(ax, (ω_grid.*2π) ./ 1e3, TOF_simple; linestyle=:solid,  label="Field-independent")

axislegend(ax, position=:rb)

save("src/WPI/test_scripts/frequencies/several_omega_$(round(rad2deg(acos(abs(μ))))).png", fig)
