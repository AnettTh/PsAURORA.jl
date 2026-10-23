using AURORA
using AURORA; c₀
using CairoMakie

## Define magnetic field
ϕ = deg2rad(240.0)
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ)
magnetic_field = DipoleMagneticField()

## Define the particle
μ = -cos(deg2rad(2))       # Almost field-aligned
L = 6.0
r0 = Cartesian(L*RE, 0.0, ϕ_eq)

E = 30e3
particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true)

## Define the plasma
λ_grid = range(deg2rad(-50), deg2rad(50), length=500)
lb_low, lb_high = 0.1, 0.5
plasma = PlasmaState(
    λ_grid,
    ϕ,
    DentonDensity(magnetic_field, L, ϕ),
    magnetic_field,
    Float64(L);
    lb_low=lb_low,
    lb_high=lb_high
    )


## Get group velocity
ω_grid=range(plasma.ω_lb[1], plasma.ω_lb[end], length=500)

## Calculate TOF
TOF_fd     = WPI_TOF(ω_grid, particle, plasma; wave_launch_time=wave_chirp, field_dependent=true)
TOF_simple = WPI_TOF(ω_grid, particle, plasma; wave_launch_time=wave_chirp, field_dependent=false)


## Make Figure
fig = Figure()
ax = Axis(fig[1,1], xlabel="ω / Ωₑ", ylabel="TOF [s]", title="Time-of-flight as a function of frequency for pitch-angle $(round(rad2deg(acos(abs(μ)))))°")

lines!(ax, ω_grid / minimum(plasma.Ω_e), TOF_fd;     linestyle=:dash, label="Field-dependent")
lines!(ax, ω_grid / minimum(plasma.Ω_e), TOF_simple; linestyle=:solid,  label="Field-independent")

axislegend(ax, position=:lt)

##
#save("src/WPI/test_scripts/frequencies/tsyg_several_omega_$(round(rad2deg(acos(abs(μ))))).png", fig)
