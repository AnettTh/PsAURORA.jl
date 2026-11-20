using AURORA
using CairoMakie
using Dates
using ProgressMeter

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

#magnetic_field = DipoleMagneticField()
magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Define the stations fieldline
ref = LTPReference(station)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Make particle
α0 = 0.5                # NOTE: If this is above the loss-cone, returns TOF=0.0
μ0 = -cos(deg2rad(α0))

E_grid = range(7e3, 40e3, length=500)

particles = Vector{ParticleState}(undef, length(E_grid))
for (i, E) in enumerate(E_grid)
    particles[i] = ParticleState(E, μ0, fieldline)
end

## Define the plasma
desnity = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma, fieldline)

##Choose frequencies to simulate for
ω_frac = range(0.1, 0.5, length=5)
ωs = ω_frac .* minimum(plasma.Ω_e)

## Calculate TOF
TOF_dependent = zeros(length(E_grid), length(ωs))  # [n_E × n_ω]
TOF_independent = zeros(length(E_grid), length(ωs))

@showprogress for (i, p) in enumerate(particles)
    TOF_dependent[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, use_boris=use_boris, field_dependent=true).tof
    TOF_independent[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, use_boris=use_boris, field_dependent=false).tof
end

## Make Figure
fig = Figure(size=(800,400))
if typeof(magnetic_field) == DipoleMagneticField
    ax = Axis(
        fig[1,1],
        xlabel="TOF [s]",
        ylabel="E [keV]",
        title="Dipole approximation with pitch-angle α = $(α0)° at $station $t",
    )
else
    ax = Axis(
        fig[1,1],
        xlabel="TOF [s]",
        ylabel="E [keV]",
        title="Tsyganenko model with pitch-angle α = $(α0)°",
    )
end

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

for (i, ω) in enumerate(ωs)
    lines!(ax, TOF_dependent[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:solid,
        label="ω / Ωₑ = $(ω_frac[i])"
    )
    lines!(ax, TOF_independent[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:dash
    )
end
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:solid,  label="field-dependent")
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:dash,   label="field-independent")

Legend(fig[1,2], ax)

##
if typeof(magnetic_field) == DipoleMagneticField
    save("src/WPI/test_scripts/field_independent_model/figures/11132026_Dipole_SeveralEnergies_α$(α0).png", fig)
else
    save("src/WPI/test_scripts/field_independent_model/figures/11132026_Tsyganenko_SeveralEnergies_α$(α0).png", fig)
end
