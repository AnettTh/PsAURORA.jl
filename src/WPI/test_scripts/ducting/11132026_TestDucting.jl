using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Define method
duct_enhancement = 0.5
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
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
ref = LTPReference(station)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Make particle
α0 = 0.5
μ0 = -cos(deg2rad(α0))

E_grid = range(7e3, 40e3, length=500)

particles = Vector{ParticleState}(undef, length(E_grid))
for (i, E) in enumerate(E_grid)
    particles[i] = ParticleState(E, μ0, fieldline)
end

## Make the plasma
ne = DentonDensity(fieldline)
ne_ducted = DentonDensity(fieldline; duct_enhancement=duct_enhancement)

plasma = PlasmaState(fieldline, ne)
plasma_ducted = PlasmaState(fieldline, ne_ducted)

## Make the wave
amplitude = ConstantAmplitude(100e-12)
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # Do not depend on density

wave = WaveState(modulation, amplitude, plasma, fieldline)
wave_ducted = WaveState(modulation, amplitude, plasma_ducted, fieldline)

## Choose frequency to simulate for
frac = 0.3
ω = frac * minimum(plasma.Ω_e)

## Calculate TOF
TOF = zeros(length(E_grid))
TOF_ducted = zeros(length(E_grid))

progress = Progress(length(particles))
Threads.@threads for i in eachindex(particles)
    TOF[i] = WPI_TOF(ω, particles[i], fieldline, plasma; wave=wave).tof
    TOF_ducted[i] = WPI_TOF(ω, particles[i], fieldline, plasma_ducted; wave=wave_ducted).tof
    next!(progress)
end

##
TOF = replace(TOF, 0.0 => NaN)
TOF_ducted = replace(TOF_ducted, 0.0 => NaN)

## Make figure
fig = Figure(size=(650, 400))
ax = Axis(
    fig[1, 1],
    xlabel="TOF [s]",
    ylabel="E [keV]",
    title="Ducted vs non-ducted time-of-flight at $t",
)

scatter!(ax, TOF, E_grid ./1e3, label="Non-ducted"; color="#E69F00")
scatter!(ax, TOF_ducted, E_grid ./1e3, label="Ducted"; color="#56B4E9")

# Add info to Legend
info = [
    "α = $(α0)°",
    "ω = $(frac)Ωₑ",
    "Enhancement: $(duct_enhancement)"
]

for i in info
    lines!(ax, [NaN], [NaN]; color=:transparent, label=i)
end

Legend(fig[1,2], ax)
save("src/WPI/test_scripts/ducting/figures/11132026_DuctVSNonDuct$(frac)Ωe.png", fig)
