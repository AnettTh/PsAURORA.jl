using AURORA
using JLD2
using Dates
using CairoMakie
using ProgressMeter

t_start = time()
stage(msg) = (println("[", lpad(round(Int, time() - t_start), 5), " s] ", msg); flush(stdout))

## Define the configuration
name = "src/WPI/data/11-13-2026/TOFCount.jld2"

git_hash = try
    readchomp(`git -C $(pkgdir(AURORA)) rev-parse --short HEAD`)
catch
    "unknown"
end

duct_enhancement = 0.0
station          = :sodankyla
WPI              = true

## Define magnetic field
stage("Defining the magnetic field")
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

## Grids to evaluate on
E_grid = range(7e3, 40e3, length=20)
α_grid = range(deg2rad(0.01), deg2rad(1.0), length=2)
μ0_grid = -cos.(α_grid)
ω_frac = range(0.1, 0.5, length=20)

## Define station fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define the particles
stage("Building the particles")

particles = Matrix{ParticleState}(undef, length(E_grid), length(μ0_grid))
for iα in eachindex(α_grid), (iE, E) in enumerate(E_grid)
    particles[iE, iα] = ParticleState(
        E,
        μ0_grid[iα],
        fieldline,
    )
end

## Define the plasma and the wave
stage("Building the plasma")
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the frequencies sampled
ωs = ω_frac .* minimum(plasma.Ω_e)

## Define the wave
a_Bw = 100e-12
amplitude = ConstantAmplitude(a_Bw)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WPI ? WaveState(modulation, amplitude, plasma, fieldline) : nothing

## Calculate TOF and FP
stage("Running the boris-mover")
TOF = fill(NaN, length(E_grid), length(μ0_grid), length(ωs))
FP = fill(NaN, 3, length(E_grid), length(μ0_grid), length(ωs))

idx = CartesianIndices(particles)
progress = Progress(length(idx))
Threads.@threads for k in eachindex(idx)
    iE, iα = Tuple(idx[k])

    result = WPI_TOF(
        ωs,
        particles[iE, iα],
        fieldline,
        plasma;
        wave=wave,
        get_footpoint=true,
    )
    TOF[iE, iα, :] = result.tof
    FP[:, iE, iα, :] .= result.footpoint

    next!(progress)
end
stage("Done: $(count(!isnan, TOF)) of $(length(TOF)) resonant")

## Save result
params = (
    # method
    WPI              = WPI,
    n                = 1,
    duct_enhancement = duct_enhancement,
    # magnetic field
    time          = string(magnetic_field.t),
    pdyn          = magnetic_field.pdyn,
    dst           = magnetic_field.dst,
    byimf         = magnetic_field.byimf,
    bzimf         = magnetic_field.bzimf,
    ps            = magnetic_field.ps,          # [rad]

    # station and traced line
    station       = station,
    r_eq          = fieldline.r_eq,             # [m]
    ϕ_eq          = fieldline.ϕ_eq,             # [rad]
    B_eq          = fieldline.B_eq,             # [T]
    r_max         = fieldline.r_max,            # [m]

    # density and plasma
    L_Denton_bin     = density.L_Denton_bin,
    lb_low           = 0.1,                     # use the values you actually passed
    lb_high          = 0.5,
    Ω_eq             = minimum(plasma.Ω_e),     # [rad/s]

    # wave
    chirp_ω0       = modulation.ω0,             # [rad/s]
    chirp_ω1       = modulation.ω1,
    chirp_duration = modulation.t_duration,     # [s]
    Bw             = amplitude.Bw,              # [T]

    # Provenance
    git_commit       = git_hash,
    date             = string(now()),
    notes            = "TOF [s]: NaN = no resonance, 0.0 = bounced. " *
                       "FP in metres, SM frame at `time` (needs SM to GEO before ENU). " *
                       "Particles built at r0 = Dipolar(L, 0, ϕ_eq); " *
                       "r_start = position_eq for every case.",
)

jldsave(name;
    TOF,                          # TOF [E, α, ω]
    FP,                           # FP [xyz, E, α, ω]
    E_grid  = collect(E_grid),
    α_grid,
    ωs,
    params,
)

@load "src/WPI/data/11-13-2026/energy-time-position.jld2"
size(TOF)       # [E, α, ω, Ψ0]
size(FP)        # [xyz, E, α, ω, Ψ0]
size(E_grid)
size(α_grid)
size(ω_frac)

ωs = ω_frac
station = params.station
t = DateTime(params.time)
##
for i in eachindex(Ψ0_grid)
    fig = Figure(size = (900, 400))

    TOF_α1 = TOF[:, 1, :, 4]
    TOF_α2 = TOF[:, 2, :, 4]

    vals = [v for s in (TOF_α1, TOF_α2) for v in s if !isnan(v) && v != 0.0]
    colorrange = (3, maximum(vals))

    Label(fig[0, 1:2], "Time-of-flight for $(station) at $t, Ψ0 = $(round(Ψ0_grid[i] / π, digits=1))π"; fontsize=18, font=:bold, tellwidth=false)

    ax1 = Axis(
        fig[1, 1],
        xlabel="ω / Ωₑ",
        ylabel="E [keV]",
        title="α = $(round(rad2deg(α_grid[1]), digits=2))°"
        )

    ax2 = Axis(
        fig[1, 2],
        xlabel="ω / Ωₑ",
        ylabel="E [keV]",
        title="α = $(round(rad2deg(α_grid[2]), digits=2))°"
    )

    hm1 = heatmap!(ax1, ωs, E_grid / 1e3, TOF_α1', colormap=:turbo, colorrange=colorrange)
    hm2 = heatmap!(ax2, ωs, E_grid / 1e3, TOF_α2', colormap=:turbo, colorrange=colorrange)

    Colorbar(fig[1, 3], hm1; label = "TOF [s]")

    save("src/WPI/test_scripts/big_data/figures/11132026-TOF/TOF$(round(Ψ0_grid[i] / π, digits=1))π.png", fig)
end
