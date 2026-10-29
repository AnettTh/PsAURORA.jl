using AURORA
using CairoMakie
using ProgressMeter

## Define a magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define a particle
μ = -cos(deg2rad(0.1))
L = 7.1
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))
E = 10e3

particle = ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true)

## Define a plasma
λ_grid = range(deg2rad(-50.0), deg2rad(50.0), length=200)
density = DentonDensity(magnetic_field, L, ϕ_eq)
plasma = PlasmaState(λ_grid, density; lb_low=0.1, lb_high=0.5)

## Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

## Find the resonance latitude
ω_start = wave_frequency(wave, 0.0)
λ_resonance = resonance_latitude(ω_start, particle, wave.plasma)[1]

## Find TOF and footpoint
time, footpoint = boris_mover(particle;
    wave=wave,
    λ_start=λ_resonance,
    n_T=100000,
    resolution=10,
    store_trajectory=false,
    count_bounces=false,
    stop_when_bouncing=true
)

##
