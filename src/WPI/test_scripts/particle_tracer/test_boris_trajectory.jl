using AURORA
using WGLMakie
using Bonito

## Define a magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
#magnetic_field = DipoleMagneticField()

# Define a particle
E = 40e3
μ = -cos(deg2rad(2))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true)

# Define a plasma
λ_grid = range(deg2rad(-50.0), deg2rad(50.0), length=200)
ne_model = DentonDensity(magnetic_field, L, ϕ_eq)

plasma = PlasmaState(λ_grid, ϕ_eq, ne_model, magnetic_field, L; lb_low=0.1, lb_high=0.5)

# Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

# Find resonance for this wave
ω_start = wave_frequency(wave, 0.0)
λ_resonance, _ = resonance_latitude(ω_start, particle, wave.plasma)

## Run Boris-mover
result = boris_mover(particle;
    wave=wave,
    λ_start=λ_resonance,
    n_T=100000,
    resolution=10,
    store_trajectory=true,
    count_bounces=true,
    stop_when_bouncing=false
)

## Make 3D figure
# Get coordinates
x, y, z = result.position[:, 1], result.position[:, 2], result.position[:, 3]
# Make the Earth sphere
u = LinRange(0, π, 50)
v = LinRange(0, 2π, 50)
x_earth = sin.(u) * cos.(v)'
y_earth = sin.(u) * sin.(v)'
z_earth = cos.(u) * ones(50)'

# Make figure
fig = Figure()
ax  = Axis3(fig[1,1],
    title  = "Particle trajectory",
    xlabel = "X",
    ylabel = "Y",
    zlabel = "Z",
    aspect=:data
)

surface!(ax, x_earth, y_earth, z_earth, color = :blue, alpha = 0.5)
lines!(ax, x[2:end]./RE, y[2:end]./RE, z[2:end]./RE)

# Display figure
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)


## Make 2D-figure
# Get coordinates — project to xz-plane
x = result.position[:, 1]
z = result.position[:, 3]

# Make figure
fig = Figure(size=(700, 600))
ax  = Axis(fig[1,1];
    title  = "Particle trajectory (xz-plane)",
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    aspect = DataAspect()
)

# Earth circle
θ_earth = range(0, 2π, length=200)
lines!(ax, cos.(θ_earth), sin.(θ_earth); color=:blue, linewidth=2, label="Earth")

# Field line
r0_field = Cartesian(Spherical(particle.L, λ_resonance, particle.ϕ_eq))
_, trace  = find_R_max(particle.magnetic_field, r0_field; store_trace=true)
lines!(ax, trace.xs./RE, trace.zs./RE; color=:gray, linestyle=:dash, label="Tsyganenko")

# Add also dipole fieldline
dipole = DipoleMagneticField()
_, trace_dip = find_R_max(dipole, r0_field; store_trace=true)
lines!(ax, trace_dip.xs./RE, trace_dip.zs./RE;
    color=:gray, alpha=0.5, linestyle=:dash, label="Dipole")

# Particle trajectory
lines!(ax, x[2:end]./RE, z[2:end]./RE; color=:red, linewidth=1, label="Trajectory")

# Start point
scatter!(ax, [x[2]/RE], [z[2]/RE]; color=:green, markersize=15, marker=:circle, label="Start")

# footpoint or endpoint
scatter!(ax, [x[end]/RE], [z[end]/RE]; color=:purple, markersize=15, marker=:circle, label="End")

axislegend(ax)

# Display with Bonito
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
