using AURORA
using AURORA; c₀
using CairoMakie
using Dates
using ProgressMeter

##
duct_enhancement = 0.0
station          = :sodankyla

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
E_grid = range(1e3, 40e3, length=100)

p_relativistic = [ParticleState(E, μ0, fieldline; relativistic=true) for E in E_grid]
p_classical = [ParticleState(E, μ0, fieldline; relativistic=false) for E in E_grid]

## Define the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
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
        color=colors[i], marker=:circle)
    scatter!(ax1, TOF_simple_relativistic[:, i], E_grid ./ 1e3;
        color=colors[i], marker=:cross)

    # Panel 2: field-dependent mode
    scatter!(ax2, TOF_matrix_classical[:, i],    E_grid ./ 1e3;
        color=colors[i], marker=:circle)
    scatter!(ax2, TOF_matrix_relativistic[:, i], E_grid ./ 1e3;
        color=colors[i], marker=:cross)
end

# Dummy lines for legend
[scatter!(ax2, [NaN], [NaN]; color=colors[i], marker=:rect, label="ω = $(ω_frac[i]) Ωₑ") for i in eachindex(ωs)]
scatter!(ax2, [NaN], [NaN]; color=:black, marker=:circle, label="Classical")
scatter!(ax2, [NaN], [NaN]; color=:black, marker=:cross,  label="Relativistic")

Legend(fig[1, 3], ax2)
xlims!(ax1, 3, 11)
xlims!(ax2, 3, 11)
##
save("src/WPI/test_scripts/relativistic/figures/11192026_RelativisticVsClassical.png", fig)
