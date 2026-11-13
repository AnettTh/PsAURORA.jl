using AURORA
using CairoMakie

## Define method
duct_enhancement = 0.0
station = :sodankyla
use_boris = false

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

# TODO: Figure out why this one doesnt run with Tsyganekofield
magnetic_field = DipoleMagneticField()
#magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Define the stations fieldline
ref = LTPReference(station)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define the particle
α0 = 2.0
μ0 = -cos(deg2rad(α0))
E_grid = [7e3, 10e3, 15e3, 20e3, 30e3, 40e3]
particles = [ParticleState(E, μ0, fieldline) for E in E_grid]

## Define the plasma
desnity = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma, fieldline)

##Choose frequencies to simulate for
ω_frac = range(0.1, 0.5, length=200)
ωs = ω_frac .* minimum(plasma.Ω_e)

## Calculate TOF
TOF_fd = [WPI_TOF(ωs, p, fieldline, plasma; wave=wave, use_boris=use_boris, field_dependent=true).tof for p in particles]
TOF_fi = [WPI_TOF(ωs, particle, fieldline, plasma; wave=wave, use_boris=use_boris, field_dependent=false).tof for p in particles]

## Make Figure
fig = Figure(size=(800, 400))
if typeof(magnetic_field) == DipoleMagneticField
    ax = Axis(
        fig[1,1],
        xlabel="ω / Ωₑ",
        ylabel="TOF [s]",
        title="Dipole approximation with pitch-angle α = $(α0)°",
    )
else
    ax = Axis(
        fig[1,1],
        xlabel="ω / Ωₑ",
        ylabel="TOF [s]",
        title="Tsyganenko model with pitch-angle α = $(α0)°",
    )
end

colors =["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

for i in eachindex(particles)
    lines!(ax, ω_frac, TOF_fd[i]; linestyle=:solid, label="$(E_grid[i] / 1e3) keV", color=colors[i])
    lines!(ax, ω_frac, TOF_fi[i]; linestyle=:dash, color=colors[i])
end
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:solid, label="field-dependent")
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:dash, label="field-independent")

Legend(fig[1,2], ax)

##
if typeof(magnetic_field) == DipoleMagneticField
    save("src/WPI/test_scripts/field_independent_model/figures/11132026_Dipole_SeveralFrequencies_α$(α0).png", fig)
else
    save("src/WPI/test_scripts/field_independent_model/figures/11132026_Tsyganenko_SeveralFrequencies_α$(α0).png", fig)
end
