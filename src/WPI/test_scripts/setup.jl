using AURORA
using Dates

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
ref = LTPReference(station)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define particle
E = 39e3
α0 = 0.5
μ0 = -cos(deg2rad(α0))
particle = ParticleState(E, μ0, fieldline)

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

Ψ0 = range(0, 2π, length=20)[2]
#waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ) for Ψ in Ψ0]
wave = WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0)

## Frequencies to evaluate
ω_frac = 0.1
ω = ω_frac * minimum(plasma.Ω_e)

#boris = [boris_mover(particle, fieldline, wave; ω=ω) for wave in waves]
boris = boris_mover(particle, fieldline, wave; ω=ω)
wpi = WPI_TOF(ω, particle, fieldline, plasma; wave=wave)


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
    ω_frac           = ω_frac,
    a_Bw             = amplitude.Bw,
    chirp_duration   = modulation.t_duration,
)

jldsave(name;
    TOF,
    FP,
    FP_enu,
    α_FP,
    E_grid=collect(E_grid),
    Ψ0_grid,
    params,
)
