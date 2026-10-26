using AURORA
using WGLMakie
using Bonito

## Define a magnetic field
ϕ = deg2rad(240.0)
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define a particle
E = 60e3
μ = -cos(deg2rad(0.01))
L = 4.0
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq)

## Define a plasma
λ_grid = range(deg2rad(-50.0), deg2rad(50.0), length=200)
ne_model = DentonDensity(magnetic_field, L, ϕ_eq)

plasma = PlasmaState(λ_grid, ϕ_eq, ne_model, magnetic_field, L; lb_low=0.1, lb_high=0.5)

## Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)


## Run Boris-mover
result = boris_mover(particle; wave=wave, n_T=1000000, resolution=5, store_trajectory=true, λ_start=-deg2rad(20))

## Get coordinates
x, y, z = result.position[:, 1], result.position[:, 2], result.position[:, 3]

## Make the Earth sphere
u = LinRange(0, π, 50)
v = LinRange(0, 2π, 50)
x_earth = sin.(u) * cos.(v)'
y_earth = sin.(u) * sin.(v)'
z_earth = cos.(u) * ones(50)'

## Make figure
fig = Figure()
ax  = Axis3(fig[1,1],
    title  = "Particle trajectory",
    xlabel = "X",
    ylabel = "Y",
    zlabel = "Z",
    aspect=:data
)

surface!(ax, x_earth, y_earth, z_earth, color = :blue, alpha = 0.5)
lines!(ax, x./RE, y./RE, z./RE)

# Display figure
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
