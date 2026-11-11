using AURORA
using JLD2
using Dates
using CairoMakie
using ProgressMeter

t_start = time()
stage(msg) = (println("[", lpad(round(Int, time() - t_start), 5), " s] ", msg); flush(stdout))

## Define the configuration
name = "src/WPI/data/11-10-2026/TOF_count.jld2"

git_hash = try
    readchomp(`git -C $(pkgdir(AURORA)) rev-parse --short HEAD`)
catch
    "unknown"
end

use_boris        = true
field_dependent  = true
relativistic     = true
WPI              = true
n                = 1
duct_enhancement = 0.0
station          = :sodankyla
wave_direction   = -1

## Define magnetic field
stage("Defining the magnetic field")
t     = DateTime("2020-12-21T02:00:00")
ref   = LTPReference(:sodankyla)
pdyn  = 2.0
dst   = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    ref;
    time=t,
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)


# Find the L-shell and equatorial position of the fieldline belonging to the station at the
# requested time
stage("Tracing the stations's fieldline")
trace_L_shell = find_L_shell(magnetic_field, ref, t)
L = trace_L_shell.L
position_eq = trace_L_shell.position_eq

## Define the particles
stage("Building the particles")
α_grid = deg2rad.([0.01, 1.0])
μ0_grid = -cos.(α_grid)

ϕ_eq = magnetic_field.ϕ
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))
E_grid = range(7e3, 40e3, length=20)

particles = Matrix{ParticleState}(undef, length(E_grid), length(μ0_grid))
@showprogress for iα in eachindex(α_grid), (iE, E) in enumerate(E_grid)
    particles[iE, iα] = ParticleState(
        E,
        μ0_grid[iα],
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=relativistic,
    )
end

## Define the plasma and the wave
stage("Building the plasma")
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
lb_low, lb_high = 0.1, 0.5
density = DentonDensity(
    magnetic_field,
    L,
    ϕ_eq;
    duct_enhancement=duct_enhancement,
)

plasma = PlasmaState(λ_grid, density; lb_low=lb_low, lb_high=lb_high)

## Define the wave
# Frequencies sampled
ωs = range(plasma.ω_lb[1], plasma.ω_lb[end], length=20)
ω_frac = range(lb_low, lb_high, length=length(ωs))

# Complete wave
a_Bw = 100e-12
amplitude = ConstantAmplitude(a_Bw)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WPI ? WaveState(
    modulation,
    amplitude,
    plasma;
    wave_direction=wave_direction
) : nothing


## Calculate TOF and FP
stage("Running the boris-mover")
TOF = fill(NaN, length(E_grid), length(μ0_grid), length(ωs))
FP  = fill(NaN, 3, length(E_grid), length(μ0_grid), length(ωs))

idx      = CartesianIndices(particles)
progress = Progress(length(idx))
Threads.@threads for k in eachindex(idx)
    iE, iα = Tuple(idx[k])

    result = WPI_TOF(
        ωs,
        particles[iE, iα],
        plasma;
        wave=wave,
        use_boris=use_boris,
        field_dependent=field_dependent,
        n=n,
        r_start=position_eq,
        get_footpoint=true,
    )
    TOF[iE, iα, :]    = result.tof
    FP[:, iE, iα, :] .= result.footpoint

    next!(progress)
end
stage("Done: $(count(!isnan, TOF)) of $(length(TOF)) resonant")

## Save result
params = (
    # method
    use_boris        = use_boris,
    relativistic     = relativistic,
    WPI              = WPI,
    field_dependent  = field_dependent,
    n                = n,
    duct_enhancement = duct_enhancement,
    # station and field
    station          = station,
    lat_deg          = ref.lat_deg,
    lon_deg          = ref.lon_deg,
    alt_m            = ref.alt_m,
    time             = string(t),
    pdyn             = pdyn,
    dst              = dst,
    byimf            = byimf,
    bzimf            = bzimf,
    ps               = magnetic_field.ps,                  # [rad]
    # Field line and plasma
    L                = L,
    ϕ_eq             = ϕ_eq,
    position_eq      = [position_eq.x, position_eq.y, position_eq.z],   # [m], SM
    lb_low           = lb_low,
    lb_high          = lb_high,
    Ω_eq             = minimum(plasma.Ω_e),                # [rad/s]
    L_Denton_bin     = density.L_Denton_bin,
    R_max            = density.R_max,                      # [m]
    # Wave
    ω_frac           = ω_frac,
    a_Bw             = a_Bw,                               # [T]
    chirp_duration   = modulation.t_duration,              # [s]
    wave_direction   = wave_direction,
    # Provenance
    git_commit       = git_hash,
    date             = string(now()),
    julia_version    = string(VERSION),
    nthreads         = Threads.nthreads(),
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
    μ0_grid,
    ωs,
    params,
)


@load "src/WPI/data/11-10-2026/TOF_count.jld2"

#size(TOF)       # [E, α, ω]
#size(FP)        # [xyz, E, α, ω]
#size(E_grid)    # [20]
#size(α_grid)    # [2]
#size(ωs)        # [20]

##
fig = Figure(size = (1100, 500))

TOF_α1 = TOF[:, 1, :]
TOF_α2 = TOF[:, 2, :]

vals = [v for s in (TOF_α1, TOF_α2) for v in s if !isnan(v) && v != 0.0]
colorrange = extrema(vals)

ax1 = Axis(
    fig[1, 1],
    xlabel="ω / Ωₑ",
    ylabel="E [keV]",
    )

ax2 = Axis(
    fig[1, 2],
    xlabel="ω / Ωₑ",
    ylabel="E [keV]",
)

hm1 = heatmap!(ax1, params.ω_frac, E_grid / 1e3, TOF_α1', colormap=:turbo)
hm2 = heatmap!(ax2, params.ω_frac, E_grid / 1e3, TOF_α2', colormap=:turbo)

Colorbar(fig[1, 3], hm1; label = "TOF [s]")
