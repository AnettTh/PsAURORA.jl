using AURORA
using AURORA; c₀

# NOTE: For several whistler mode wave angular frequencies
## Define the particle
μ = -0.99       # Almost field-aligned
L = 6.0
r0 = [L*RE, 0.0, 0.0]

E = 3e4
particle = ParticleState(E, μ, r0, dipole_field)

# TODO: λ_grid here should be found by using L--integrate this into PlasmaParameters!
## Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
n_e0 = 1.8e7        # 18/cc, from Hsieh 2022

plasma = PlasmaParameters(λ_grid, n_e0, dipole_field, L)


## Get parallel velocity
v_parallel = parallel_velocity(particle, λ_grid)


## Get group velocity
ω_grid=range(1e3, 4.74e3, length=500)
v_g_grid = group_velocity_whistler_wave(ω_grid)


## Find the field-dependent time-of-flight
TOF_fd = WPI_TOF_field_dependent(ω_grid, plasma, particle, t0)


## Find the field-independent time-of-flight
v_parallel_lc = particle.v * cos(particle.α_lc)
TOF = [WPI_TOF(deg2rad(20), v_g, v_parallel_lc, particle.L*RE, 0, particle.E_eV, particle.α_lc) for v_g in v_g_grid]


## Make Figure
fig = Figure()
ax = Axis(fig[1,1], xlabel="ω [kHz]", ylabel="time-of-flight [s]")
l1 = lines!(ax, ω_grid./1e3, TOF_fd, label="Chen")
l2 = lines!(ax, ω_grid./1e3, TOF, label="Miyoshi-Saito")

fig[1, 2] = Legend(fig, ax, "Method used")

fig
