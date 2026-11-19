using AURORA
using CairoMakie
using ProgressMeter
using Dates

##
duct_enhancement = 0.0
station = :sodankyla

## Field, station line, particle and plasma
t = DateTime("2020-12-21T02:00:00")
magnetic_field = TsyganenkoMagneticField(t; pdyn=2.0, dst=-28.0, byimf=0.0, bzimf=2.0)
fieldline = StationFieldLine(magnetic_field, LTPReference(station), t)

E_grid = range(7e3, 40e3, length=300)
α0 = 0.5
μ0 = -cosd(α0)
particles = [ParticleState(E, μ0, fieldline) for E in E_grid]

density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
amplitude = ConstantAmplitude(100e-12)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WaveState(modulation, amplitude, plasma, fieldline)

## Choose frequencies to run for
ω_frac = range(0.1, 0.5, length=50)
ωs = ω_frac .* minimum(plasma.Ω_e)

## Calculate TOF
TOF = zeros(length(E_grid), length(ωs))
progress = Progress(length(particles))

Threads.@threads :greedy for i in eachindex(particles)
    TOF[i, :] = WPI_TOF(ωs, particles[i], fieldline, plasma; wave=wave).tof
    next!(progress)
end


## Make Figure
fig = Figure(size=(700,500))
ax = Axis(
    fig[1,1],
    xlabel="ω/Ωₑ",
    ylabel="E [keV]",
    title="Time-of-flight for $station station at $t   (α₀ = $(α0)°)"
)

hm = heatmap!(ax, ω_frac, E_grid ./ 1e3, TOF'; colormap=:turbo, colorrange=(2.5,4.3))
Colorbar(fig[1,2], hm; label="TOF [s]")
##
save("src/WPI/test_scripts/TOF/figures/11192026_EnergyHeatmapTOF.png", fig)
