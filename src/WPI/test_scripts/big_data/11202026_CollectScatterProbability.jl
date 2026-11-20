using AURORA
using JLD2
using ProgressMeter
using Dates

## Configurations; rest is default
duct_enhancement = 0.0
station = :sodankyla
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0
chirp_duration = 2.0
a_Bw = 100e-12


## Grids with TOF(E, μ0, ω, Ψ0) = (6, 250, 3, 48)
E_grid = [7e3, 10e3, 15e3, 20e3, 30e3, 40e3]
α0_grid = range(0.1, 10, length=250)
μ0_grid = -cosd.(α0_grid)
ω_frac = [0.1, 0.38, 0.5]
Ψ0_grid = range(0, 2π, length=48)[1:end-1]

## Output
TOF = fill(NaN, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))
FP = fill(NaN, 3, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))
outfile = "src/WPI/data/11-20-2026/CollectScatterProbability2.jld2"

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
        α0_grid,
        ω_frac,
        Ψ0_grid,
        params,
    )
end

## Define magnetic field
magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Make the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the frequencies sampeled
ωs = ω_frac .* minimum(plasma.Ω_e)

## Make the wave modulation for the full lower-band whistler mode chorus wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], chirp_duration)
amplitude = ConstantAmplitude(a_Bw)

## Store parameters of interest
params = (
    # provenance
    git_commit   = git_hash,
    git_dirty    = git_dirty,
    date_started = date_started,

    # magnetic field
    time  = string(magnetic_field.t),
    pdyn  = magnetic_field.pdyn,
    dst   = magnetic_field.dst,
    byimf = magnetic_field.byimf,
    bzimf = magnetic_field.bzimf,
    ps    = magnetic_field.ps,          # [rad]

    # station and traced line
    station = station,
    r_eq    = fieldline.r_eq,             # [m]
    ϕ_eq    = fieldline.ϕ_eq,             # [rad]
    B_eq    = fieldline.B_eq,             # [T]
    r_max   = fieldline.r_max,            # [m]

    # density and plasma
    duct_enhancement = duct_enhancement,
    L_Denton_bin     = density.L_Denton_bin,
    lb_low           = plasma.ω_lb[1],
    lb_high          = plasma.ω_lb[end],
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
    particles[iE, iμ] = ParticleState(E_grid[iE], μ0_grid[iμ], fieldline)
end

## Running the simulation for each wave- and particle combination
for (iΨ, Ψ0) in enumerate(Ψ0_grid)
    # Make wave with current Ψ0
    wave = WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0)

    # Track the progress
    progress_run = Progress(length(index_pairs); desc="Ψ₀ = $(round(rad2deg(Ψ0), digits=0))°: ")

    # Run the TOF calculation for each of the energy- and pitch-angle pairs
    Threads.@threads :greedy for k in eachindex(index_pairs)
        iE, iμ = Tuple(index_pairs[k])
        p = particles[iE, iμ]
        result = WPI_TOF(ωs, p, fieldline, plasma; wave=wave)

        TOF[iE, iμ, :, iΨ] = result.tof
        FP[:, iE, iμ, :, iΨ] = result.footpoint

        next!(progress_run)
    end
end

save_results(outfile, TOF, FP, params)
