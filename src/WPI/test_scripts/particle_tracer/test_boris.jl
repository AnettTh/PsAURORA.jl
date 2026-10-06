using AURORA

## Define a magnetic field
ϕ = deg2rad(120.0)
ϕ_eq = deg2rad(120.0)
magnetic_field = TsyganenkoMagneticField(ϕ)

## Define a particle
E = 30e3
μ = cos(deg2rad(2))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ))      # TODO: Make something that converts from spherical to cartesian

particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq)


## Run Boris-mover
result = boris_mover(particle; n_T=100000, store_trajectory=false)
