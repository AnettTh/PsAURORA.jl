using AURORA
using Dates
using ProgressMeter
using CairoMakie
using JLD2

## Configuration
duct_enhancement = 0.0
station = :sodankyla
name = "src/WPI/data/11-26-2026/FootpointForFrequencies.jld2"

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

Ψ0_grid = range(0.0, 2π, length=25)[1:end-1]
waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Frequencies to evaluate
ω_frac = range(0.1, 0.5, length=100)
ω = ω_frac * minimum(plasma.Ω_e)

nE, nω, nΨ = length(E), length(ω), length(Ψ0_grid)
TOF = fill(NaN, nE, nω, nΨ)
α_FP = fill(NaN, nE, nω, nΨ)
FP = fill(NaN, 3, nE, nω, nΨ)

# Run every pair of energies and phases
pairs_EΨ = CartesianIndices((nE, nΨ))
progress = Progress(length(pairs_EΨ))
Threads.@threads :greedy for k in eachindex(pairs_EΨ)
    iE, iΨ = Tuple(pairs_EΨ[k])
    result = WPI_TOF(collect(ω), particle[iE], fieldline, plasma; wave=waves[iΨ])
    if !isnothing(result.footpoint)
        TOF[iE, :, iΨ] = result.tof
        α_FP[iE, :, iΨ] = result.α_footpoint
        FP[:, iE, :, iΨ] .= result.footpoint
    end
    next!(progress)
end

FP_enu = fill(NaN, size(FP))
for iΨ in 1:nΨ, iω in 1:nω, iE in 1:nE
    any(isnan, @view FP[:, iE, iω, iΨ]) && continue
    FP_enu[:, iE, iω, iΨ] .= FP_to_ENU(FP[:, iE, iω, iΨ], ref, t)
end

params = (
    duct_enhancement = duct_enhancement,
    station          = station,
    α0               = α0,
    time             = string(t),
    pdyn             = pdyn,
    dst              = dst,
    byimf            = byimf,
    bzimf            = bzimf,
    ps               = magnetic_field.ps,
    Ω_eq             = minimum(plasma.Ω_e),
    L_Denton_bin     = density.L_Denton_bin,
    R_max            = density.R_max,
    a_Bw             = amplitude.Bw,
    chirp_duration   = modulation.t_duration,
)

jldsave(name;
    TOF,
    FP,
    FP_enu,
    α_FP,
    E_grid=collect(E),
    ω_frac,
    Ψ0_grid,
    params,
)
