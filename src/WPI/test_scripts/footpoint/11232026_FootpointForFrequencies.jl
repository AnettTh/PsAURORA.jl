using AURORA
using Dates
using ProgressMeter
using CairoMakie

## Configuration
duct_enhancement = 0.0
station = :sodankyla
name = "src/WPI/data/11232026/FootpointForFrequencies.jl"

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
E = [7e3, 10e3, 15e3, 20e3, 30e3, 40e3]
α0 = 1.5
μ0 = -cos(deg2rad(α0))

particle = [ParticleState(e, μ0, fieldline) for e in E]

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

Ψ0_grid = range(0.0, 2π, length=13)[1:end-1]
waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Frequencies to evaluate
ω_frac = range(0.1, 0.5, length=10)
ω = ω_frac * minimum(plasma.Ω_e)

nE, nΨ = length(E), length(Ψ0_grid)
α_FP = fill(NaN, nE, nΨ)
FP = fill(NaN, 3, nE, nΨ)

# Run every pair of energies and phases
pairs_EΨ = CartesianIndices((nE, nΨ))
progress = Progress(length(pairs_EΨ))
Threads.@threads :greedy for k in eachindex(pairs_EΨ)
    iE, iΨ = Tuple(pairs_EΨ[k])
    result = boris_mover(particle[iE], fieldline, waves[iΨ]; ω=ω)
    if !isnothing(result.footpoint)
        α_FP[iE, iΨ] = result.α_footpoint
        FP[:, iE, iΨ] .= result.footpoint
    end
    next!(progress)
end

fp_enu = fill(NaN, size(FP))
for iΨ in 1:nΨ, iE in 1:nE
    any(isnan, @view FP[:, iE, iΨ]) && continue
    fp_enu[:, iE, iΨ] .= FP_to_ENU(FP[:, iE, iΨ], ref, t)
end
east_km = fp_enu[1, :, :] ./ 1e3

E_keV = collect(E) ./ 1e3
E_all = repeat(E_keV, 1, nΨ)
valid = isfinite.(α_FP)
lim = maximum(east_km[valid])
