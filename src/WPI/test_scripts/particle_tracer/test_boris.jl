using AURORA
using WGLMakie
using Bonito

## Define a magnetic field
ϕ = deg2rad(240.0)
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define a particle
E = 30e3
μ = -cos(deg2rad(2))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))      # TODO: Make something that converts from spherical to cartesian

particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq; precomputed_magnetic_field=true)


## Run Boris-mover
result = boris_mover(particle; n_T=100000, store_trajectory=true, λ_start=-deg2rad(20))

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

## Display figure
app = Bonito.App(fig)
#server = Bonito.Server(app, "0.0.0.0", 8888)
