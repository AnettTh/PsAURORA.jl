using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Define model
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T03:00:00")
ref = LTPReference(:sodankyla)
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define the particle
α_deg = 0.5
μ0 = -cos(deg2rad(α_deg))
E_grid = range(20e3, 30e3, length=100)

particles = [ParticleState(E, μ0, fieldline) for E in E_grid]


## Define the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)
wave = WaveState(modulation, amplitude, plasma, fieldline)

ω_frac = range(0.1, 0.5, length=5)
ωs = ω_frac * minimum(plasma.Ω_e)

##
TOF_n1 = zeros(length(E_grid), length(ωs))
TOF_n2 = zeros(length(E_grid), length(ωs))

@showprogress for (i, p) in enumerate(particles)
    TOF_n1[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, n=1, use_boris=true).tof
    TOF_n2[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, n=2, use_boris=true).tof
end


##
fig = Figure(size=(1200, 500))

Label(fig[0, 1:3],
    "Wave-particle interaction time-of-flight";
    fontsize=16, font=:bold
)

ax1 = Axis(fig[1, 1];
    xlabel = "Time [s]",
    ylabel = "E [keV]",
    title  = "n = 1"
)

ax2 = Axis(fig[1, 2];
    xlabel = "Time [s]",
    ylabel = "E [keV]",
    title  = "n = 2"
)

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

for (i, ω) in enumerate(ωs)
    tof1 = replace(TOF_n1[:, i], 0.0 => NaN)
    tof2 = replace(TOF_n2[:, i], 0.0 => NaN)
    scatter!(ax1, tof1, E_grid./1e3, label="ω = $(ω_frac[i]) Ωₑ"; color=colors[i])
    scatter!(ax2, tof2, E_grid./1e3; color=colors[i])
end

linkyaxes!(ax1, ax2)
linkxaxes!(ax1, ax2)

Legend(fig[1, 3], ax1)
save("src/WPI/test_scripts/resonance/figures/11192026_CompareTOD_n1n2.png", fig)
