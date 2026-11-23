using AURORA
using AURORA; c₀
using CairoMakie
using Dates
using ProgressMeter
using JLD2

##
duct_enhancement = 0.0
station = :sodankyla
name = "11-23-2026/RelativisticCorrection.JLD2"

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

## Station fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define the particles
α0 = 0.5
μ0 = -cos(deg2rad(α0))
E_grid = range(7e3, 40e3, length=1000)

p_relativistic = [ParticleState(E, μ0, fieldline; relativistic=true) for E in E_grid]
p_classical = [ParticleState(E, μ0, fieldline; relativistic=false) for E in E_grid]

## Define the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

wave = WaveState(modulation, amplitude, plasma, fieldline)

## Choose frequencies
ω_frac = range(0.1, 0.5, length=5)
ωs = minimum(plasma.Ω_e) * ω_frac

## Calculate TOF
TOF_matrix_relativistic = zeros(length(E_grid), length(ωs))
TOF_simple_relativistic = zeros(length(E_grid), length(ωs))

@showprogress for (i, p) in enumerate(p_relativistic)
    TOF_matrix_relativistic[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, field_dependent=true).tof
    TOF_simple_relativistic[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, field_dependent=false).tof
end

##
TOF_matrix_classical = zeros(length(E_grid), length(ωs))
TOF_simple_classical = zeros(length(E_grid), length(ωs))

@showprogress for (i, p) in enumerate(p_classical)
    TOF_matrix_classical[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, field_dependent=true).tof
    TOF_simple_classical[i, :] = WPI_TOF(ωs, p, fieldline, plasma; wave=wave, field_dependent=false).tof
end

params = (
    # magnetic field
    t       = string(magnetic_field.t),
    pdyn    = magnetic_field.pdyn,
    dst     = magnetic_field.dst,
    byimf   = magnetic_field.byimf,
    bzimf   = magnetic_field.bzimf,
    ps      = magnetic_field.ps,            # [rad]

    # station and traced line
    station = station,
    r_eq    = fieldline.r_eq,               # [m]
    ϕ_eq    = fieldline.ϕ_eq,               # [rad]
    B_eq    = fieldline.B_eq,               # [T]
    r_max   = fieldline.r_max,              # [m]

    # density and plasma
    duct_enhancement = duct_enhancement,
    L_Denton_bin     = density.L_Denton_bin,
    Ω_eq             = minimum(plasma.Ω_e), # [rad/s]

    chirp_ω0       = modulation.ω0,         # [rad/s]
    chirp_ω1       = modulation.ω1,
    chirp_duration = modulation.t_duration, # [s]
    Bw             = amplitude.Bw,          # [T]
    Ψ0             = wave.Ψ0,               # [rad]

    α0 = α0,

    notes = "TOF [s]: NaN = no resonance, 0.0 = bounced.",
)

jldsave(name;
    TOF,
    E_grid = collect(E_grid),               # [eV]
    ω_frac = collect(ω_frac),
    params,
)

##
fig = Figure(size=(1100, 400))
Label(fig[0, 1:2],
    "Energy as a function of `time-of-flight` for $(station) at $t, α₀ = $(α0)°";
    fontsize=16,
    font=:bold,
    tellwidth=false,
    halign=:center)

# Panel 1: Simple mode (field-independent), both classical and relativistic
ax1 = Axis(fig[1, 1],
    xlabel = "time-of-flight [s]",
    ylabel = "E [keV]",
    title  = "Field-independent"
)

# Panel 2: Field-dependent mode, both classical and relativistic
ax2 = Axis(fig[1, 2],
    xlabel = "time-of-flight [s]",
    ylabel = "E [keV]",
    title  = "Field-dependent"
)

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

for (i, ω) in enumerate(ωs)
    # Panel 1: simple mode
    scatter!(ax1, TOF_simple_classical[:, i], E_grid ./ 1e3;
        color=colors[i])
    scatter!(ax1, TOF_simple_relativistic[:, i], E_grid ./ 1e3;
        color=colors[i], alpha=0.5)

    # Panel 2: field-dependent mode
    scatter!(ax2, TOF_matrix_classical[:, i],    E_grid ./ 1e3;
        color=colors[i])
    scatter!(ax2, TOF_matrix_relativistic[:, i], E_grid ./ 1e3;
        color=colors[i], alpha=0.5)
end

# Dummy lines for legend
[scatter!(ax2, [NaN], [NaN]; color=colors[i], marker=:rect, label="ω = $(ω_frac[i]) Ωₑ") for i in eachindex(ωs)]
scatter!(ax2, [NaN], [NaN]; color=:black, marker=:circle, label="Classical")
scatter!(ax2, [NaN], [NaN]; color=:black, marker=:cross,  label="Relativistic")

Legend(fig[1, 3], ax2)
xlims!(ax1, 3, 11)
xlims!(ax2, 3, 11)
##
save("src/WPI/test_scripts/relativistic/figures/11232026_RelativisticVsClassical.png", fig)
