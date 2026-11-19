using AURORA
using Dates
using ProgressMeter

## Configuration
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define particle
E = 30e3
α0_grid = range(0.1, 10, length=10)
μ0_grid = -cosd.(α0_grid)
particles = [ParticleState(E, μ0, fieldline) for μ0 in μ0_grid]

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

Ψ0_grid = range(0, 2π, 11)[1:end-1]
waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Frequencies to evaluate
ω_frac = 0.3
ω =  ω_frac * minimum(plasma.Ω_e)

## Get the results
α_FP = Matrix{Any}(nothing, length(μ0_grid), length(Ψ0_grid))
@showprogress for (i, p) in enumerate(particles), (j, w) in enumerate(waves)
    α_FP[i, j] = boris_mover(p; wave=w).α_footpoint
end
