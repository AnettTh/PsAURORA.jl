using AURORA
using JLD2
using ProgressMeter
using Dates


## Configurations; rest is default
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Grids with TOF(E, μ0, ω, Ψ0) = (150, 4, 150, 3)
E_grid = range(7e3, 40e3, length=150)
α_grid = [0.5, 1.5, 5.0]
μ0_grid = -cos.(deg2rad.(α_grid))
ω_frac = range(0.1, 0.5, length=150)
Ψ0_grid = range(0, 2π, length=5)[1:end-1]

## Output
TOF = fill(NaN, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))
FP = fill(NaN, 3, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))

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


## Make all combinations of energies and pitch-angles
index_pairs = CartesianIndices((length(E_grid), length(μ0_grid)))

## Build the grid of particles for each combination of E and μ
particles = Matrix{Union{Nothing, ParticleState}}(nothing, size(index_pairs)...)
progress  = Progress(length(index_pairs); desc="Building particles: ")

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
    progress_run = Progress(length(index_pairs); desc="Ψ₀ = $(round(rad2deg(Ψ0), digits=0))°: ")

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

        next!(progress_run)
    end
end

##
ref = LTPReference(station)
fp_enu = fill(NaN, size(FP))

for iΨ in Base.axes(FP, 5), iω in Base.axes(FP, 4), iμ in Base.axes(FP, 3), iE in Base.axes(FP, 2)
    fp_enu[:, iE, iμ, iω, iΨ] .= FP_to_ENU(FP[:, iE, iμ, iω, iΨ], ref, t)
end

E_keV = collect(E_grid) ./ 1e3

east  = fp_enu[1, :, :, :, :] ./ 1e3      # n_E × n_μ × n_ω × n_Ψ, in km
north = fp_enu[2, :, :, :, :] ./ 1e3

## Save the results
outfile = "src/WPI/data/11-26-2026/energy-time-position.jld2"
mkpath(dirname(outfile))                          # JLD2 does not create folders

params = (
    # station and field
    station = station,
    time    = string(t),
    pdyn    = pdyn,
    dst     = dst,
    byimf   = byimf,
    bzimf   = bzimf,
    ps      = magnetic_field.ps,                  # [rad]

    # traced line
    r_eq  = fieldline.r_eq,                       # [m]
    ϕ_eq  = fieldline.ϕ_eq,                       # [rad]
    B_eq  = fieldline.B_eq,                       # [T]
    r_max = fieldline.r_max,                      # [m]

    # density and plasma
    duct_enhancement = duct_enhancement,
    L_Denton_bin     = density.L_Denton_bin,
    R_max_density    = density.R_max,             # [m]
    Ω_eq             = minimum(plasma.Ω_e),       # [rad/s]

    # wave
    chirp_ω0       = modulation.ω0,               # [rad/s]
    chirp_ω1       = modulation.ω1,
    chirp_duration = modulation.t_duration,       # [s]
    Bw             = amplitude.Bw,                # [T]

    axes  = "TOF: [E, α0, ω, Ψ0]; FP, fp_enu: [xyz or ENU, E, α0, ω, Ψ0]",
    notes = "TOF [s]: NaN = no resonance, 0.0 = bounced. FP in m, SM frame at `time`. " *
            "fp_enu in m relative to the station. α0 in degrees at the magnetic equator. " *
            "Only n = 1 interacts in the mover.",
)

jldsave(outfile;
    TOF,
    FP,
    fp_enu,
    E_grid  = collect(E_grid),                    # [eV]
    α_grid  = collect(α_grid),                    # [°]
    ω_frac  = collect(ω_frac),                    # fraction of Ω_eq
    ωs      = collect(ωs),                        # [rad/s]
    Ψ0_grid = collect(Ψ0_grid),                   # [rad]
    params,
)
