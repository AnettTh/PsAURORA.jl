using AURORA
using CairoMakie

## Define a magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define a particle
E = 30e3
μ = -cos(deg2rad(2))
L = 6.5
r0 = Cartesian(L*RE, 0.0, ϕ_eq)

particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq)

## Define a plasma
λ_grid = range(deg2rad(-20.0), deg2rad(20.0), length=100)
ne_model = DentonDensity(magnetic_field, L, ϕ_eq)

plasma = PlasmaState(λ_grid, ϕ_eq, ne_model, magnetic_field, L; lb_low=0.1, lb_high=0.5)


## Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)


## Make a wave field
time = 0.0:1:10     # Current time
s = 1.0             # Current position
Ψ = 0.0             # Current wave phase

B0 = particle.B0    # Initial position magnetic field
λ = deg2rad(0.0)    # Latitude


##
for t in time
    δB, δE = wave_fields(wave, t, s, Ψ, B0, λ)
    @show δB, δE
end
