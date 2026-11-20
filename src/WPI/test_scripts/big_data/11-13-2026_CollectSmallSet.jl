using AURORA
using JLD2
using ProgressMeter
using Dates
using CairoMakie

## Configurations; rest is default
duct_enhancement = 0.0
station = :sodankyla

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

## Grids with TOF(E, μ0, ω, Ψ0) = (10, 8, 5, 4)
E_grid = range(25e3, 30e3, length=10)
α_grid = range(deg2rad(0.8), deg2rad(1.0), length=8)
μ0_grid = -cos.(α_grid)
ω_frac = range(0.3, 0.35, length=5)
Ψ0_grid = collect(0:π/2:7π/4)

## Output
TOF = fill(NaN, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))
FP = fill(NaN, 3, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))
outfile = "src/WPI/data/11-13-2026/SmallTestSet.jld2"

mkpath(dirname(outfile))

git_hash = try
    readchomp(`git -C $(pkgdir(AURORA)) rev-parse --short HEAD`)
catch
    "unknown"
end
git_dirty = try
    !isempty(readchomp(`git -C $(pkgdir(AURORA)) status --porcelain`))
catch
    missing
end
date_started = string(now())


function save_results(file, TOF, FP, params)
    jldsave(
        file;
        TOF,
        FP,
        E_grid,
        α_grid,
        ω_frac,
        Ψ0_grid,
        params,
    )
end

## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Make the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the frequencies sampeled
ωs = ω_frac .* minimum(plasma.Ω_e)

## Make the wave modulation for the full lower-band whistler mode chorus wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

## Store parameters of interest
params = (
    # provenance
    git_commit    = git_hash,
    git_dirty     = git_dirty,
    date_started  = date_started,

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
    duct_enhancement = duct_enhancement,
    L_Denton_bin     = density.L_Denton_bin,
    lb_low           = 0.1,                     # use the values you actually passed
    lb_high          = 0.5,
    Ω_eq             = minimum(plasma.Ω_e),     # [rad/s]

    # wave
    chirp_ω0       = modulation.ω0,             # [rad/s]
    chirp_ω1       = modulation.ω1,
    chirp_duration = modulation.t_duration,     # [s]
    Bw             = amplitude.Bw,              # [T]

    # plus the adaptive-step settings: coarse resolution and the B threshold
    notes = "TOF [s]: NaN = no resonance, 0.0 = bounced. Positions in SM at `time`. " *
            "Only n = 1 interacts in the mover.",
)

## Make all combinations of energies and pitch-angles
index_pairs = CartesianIndices((length(E_grid), length(μ0_grid)))

## Build the grid of particles for each combination of E and μ
particles = Matrix{Union{Nothing, ParticleState}}(nothing, size(index_pairs)...)
Threads.@threads :greedy for k in eachindex(index_pairs)
        iE, iμ = Tuple(index_pairs[k])
        particles[iE, iμ] = ParticleState(
            E_grid[iE],
            μ0_grid[iμ],
            fieldline;
        )
end

## Running the simulation for each wave- and particle combination
for (iΨ, Ψ0) in enumerate(Ψ0_grid)
    # Make wave with current Ψ0
    wave = WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0)

    # Track the progress
    progress = Progress(length(index_pairs); desc="Ψ₀ = $(round(rad2deg(Ψ0), digits=0))°: ")

    # Run the TOF calculation for each of the energy- and pitch-angle pairs
    Threads.@threads :greedy for k in eachindex(index_pairs)
        iE, iμ = Tuple(index_pairs[k])
        p = particles[iE, iμ]
        result = WPI_TOF(
            ωs,
            p,
            fieldline,
            plasma;
            wave=wave,
            get_footpoint=true,
        )

        TOF[iE, iμ, :, iΨ] = result.tof
        FP[:, iE, iμ, :, iΨ] = result.footpoint

        next!(progress)
    end
end

save_results(outfile, TOF, FP, params)

## Inspect data
data = load("src/WPI/data/11-13-2026/SmallTestSet.jld2")

E_grid = data["E_grid"] ./ 1e3
α_grid = rad2deg.(data["α_grid"])
ω_frac = data["ω_frac"]
Ψ0_grid = data["Ψ0_grid"]
TOF = data["TOF"]

## Make one figure for each pitch-angle
for iα in eachindex(α_grid)
    # Make figure for energy
    fig = Figure()
    Label(fig[0, 1:3], "Time-of-flight for α = $(@sprintf("%.3f", α_grid[iα]))°";
        fontsize  = 18,
        font      = :bold,
        tellwidth = false,
    )

    # One for each Ψ₀
    ax1 = Axis(fig[1,1], ylabel="E [keV]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[1]), digits=0))°")
    ax2 = Axis(fig[1,2], title="Ψ₀=$(round(rad2deg(Ψ0_grid[2]), digits=0))°")
    ax3 = Axis(fig[2,1], xlabel="ω/Ωₑ", ylabel="E [keV]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[3]), digits=0))°")
    ax4 = Axis(fig[2,2], xlabel="ω/Ωₑ", title="Ψ₀=$(round(rad2deg(Ψ0_grid[4]), digits=0))°")

    # Same scale for the colorbar
    vals = [v for k in 1:4 for v in TOF[:, 4, :, k] if !isnan(v) && v != 0.0]
    colorrange = (minimum(vals) - 0.1, maximum(vals))

    # Make heatmap
    hm1 = heatmap!(ax1, ω_frac, E_grid, TOF[:, iα, :, 1]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm2 = heatmap!(ax2, ω_frac, E_grid, TOF[:, iα, :, 2]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm3 = heatmap!(ax3, ω_frac, E_grid, TOF[:, iα, :, 3]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm4 = heatmap!(ax4, ω_frac, E_grid, TOF[:, iα, :, 4]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)

    # Add the colorbar
    Colorbar(fig[1:2, 3], hm1; label="TOF [s]")

    save("src/WPI/test_scripts/big_data/figures/11132026-energy/TOF_$(round(α_grid[iα], digits=3)).png", fig)
end

## Make one figure for each energy
for iE in eachindex(E_grid)
    ## Make figure for pitch-angle
    fig = Figure()
    Label(fig[0, 1:3], "Time-of-flight for E = $(@sprintf("%.3f", E_grid[iα])) keV";
        fontsize  = 18,
        font      = :bold,
        tellwidth = false,
    )
    # One for each Ψ₀
    ax1 = Axis(fig[1,1], ylabel="α [°]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[1]), digits=0))°")
    ax2 = Axis(fig[1,2], title="Ψ₀=$(round(rad2deg(Ψ0_grid[2]), digits=0))°")
    ax3 = Axis(fig[2,1], xlabel="ω/Ωₑ", ylabel="α [°]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[3]), digits=0))°")
    ax4 = Axis(fig[2,2], xlabel="ω/Ωₑ", title="Ψ₀=$(round(rad2deg(Ψ0_grid[4]), digits=0))°")

    # Same scale for the colorbar
    vals = [v for k in 1:4 for v in TOF[5, :, :, k] if !isnan(v) && v != 0.0]
    colorrange = (minimum(vals) - 0.1, maximum(vals))

    # Make heatmap
    hm1 = heatmap!(ax1, ω_frac, α_grid, TOF[iE, :, :, 1]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm2 = heatmap!(ax2, ω_frac, α_grid, TOF[iE, :, :, 2]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm3 = heatmap!(ax3, ω_frac, α_grid, TOF[iE, :, :, 3]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm4 = heatmap!(ax4, ω_frac, α_grid, TOF[iE, :, :, 4]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)

    # Add the colorbar
    Colorbar(fig[1:2, 3], hm1; label="TOF [s]")

    save("src/WPI/test_scripts/big_data/figures/11132026-pitchangle/TOF_$(round(E_grid[iE], digits=1))keV.png", fig)
end
