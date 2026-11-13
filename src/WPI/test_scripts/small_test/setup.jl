using AURORA
using Dates

duct_enhancement = 0.0
station          = :sodankyla

## Define magnetic field
t     = DateTime("2020-12-21T02:00:00")
pdyn  = 2.0
dst   = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    t;
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)

## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define particle
E0 = 30e3
μ0 = -cos(deg2rad(1.0))
particle = ParticleState(E0, μ0, fieldline; relativistic=true)

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

wave = WaveState(modulation, amplitude, plasma, fieldline)


## Run the TOF calculations
ω = 0.3 * minimum(plasma.Ω_e)
WPI_TOF(ω, particle, fieldline, plasma; wave=wave, use_boris=true)
