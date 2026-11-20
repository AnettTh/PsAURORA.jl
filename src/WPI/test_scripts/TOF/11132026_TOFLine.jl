using AURORA
using CairoMakie
using ProgressMeter
using JLD2
using Dates

## Define method
name = "src/WPI/data/11-19-2026/TOFLine_Ducted0.5.jld2"
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

## Define the particles
α0 = 0.5
μ0 = -cos(deg2rad(α0))
E_grid = range(7e3, 40e3, length=1000)

particles = Vector{ParticleState}(undef, length(E_grid))
for (i, E) in enumerate(E_grid)
    particles[i] = ParticleState(E, μ0, fieldline)
end

## Define the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
amplitude = ConstantAmplitude(100e-12)
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WaveState(modulation, amplitude, plasma, fieldline)

## Choose frequencies to run for
ω_frac = range(0.1, 0.5, length=5)
ωs = ω_frac .* minimum(plasma.Ω_e)

## Calculate TOF
TOF = zeros(length(E_grid), length(ωs))

progress = Progress(length(particles))
Threads.@threads for i in eachindex(particles)
    TOF[i, :] = WPI_TOF(ωs, particles[i], fieldline, plasma; wave=wave).tof
    next!(progress)
end


## Save the data
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

## Make the figure
fig = Figure(size=(900, 500))
ax = Axis(
    fig[1,1],
    xlabel="TOF [s]",
    ylabel="E [keV]",
    title="Time-of-flight with α = $(α0)°",
)

colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]


for (i, ω) in enumerate(ωs)
    tof = replace(TOF[:, i], 0.0 => NaN)
    scatter!(ax, tof, E_grid ./ 1e3; color = colors[i], label="ω/Ωₑ=$(ω_frac[i])")
end

information = [
    "$(t)",
    "pdyn = $(pdyn)",
    "dst = $(dst)",
    "byimf = $(byimf)",
    "bzimf = $(bzimf)",
    "$(duct_enhancement) duct enhancement",
]
for i in eachindex(information)
    lines!(ax, [NaN], [NaN], label=information[i], color=:transparent)
end

Legend(fig[1,2], ax)


#save("src/WPI/test_scripts/TOF/figures/11192026_TOFLine_Ducted$(duct_enhancement).png", fig)
