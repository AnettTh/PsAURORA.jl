using AURORA

## Define a magnetic field
ϕ = deg2rad(120)
magnetic_field = TsyganenkoMagneticField(ϕ)

## Define a particle
E = 30e3
μ = -cos(deg2rad(2))
L = 6.5
r0 = Cartesian(L*RE, 0.0, 0.0)      # TODO: Make something that converts from spherical to cartesian

particle = ParticleState(E, μ, r0, magnetic_field)

## Define a plasma
λ_grid = range(deg2rad(0.0), deg2rad(20.0), length=100)
ne_model = DentonDensity(magnetic_field, L, ϕ)

plasma = PlasmaState(λ_grid, ϕ, ne_model, magnetic_field, L)


## Define a wave
θ = deg2rad(0.0)
t_duration = 1.0
A_Bw = 1e-12
wave = DemekhovWave(plasma, t_duration, θ, A_Bw)
time = range(0.0, t_duration, length=length(plasma.ω_lb))


##
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time [s]", ylabel="Frequency [rad/s]")

lines!(ax, time, [wave_frequency(wave, t) for t in time])
fig
