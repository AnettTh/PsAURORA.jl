using AURORA
using CairoMakie
using JLD2
using Dates
using ProgressMeter

## Define configuration
name = "src/WPI/data/11-09-2026/footpoint_and_phase_ducted.jld2"

git_hash = try
    readchomp(`git -C $(pkgdir(AURORA)) rev-parse --short HEAD`)
catch
    "unknown"
end

use_boris        = true             # Use particle tracer or direct integration
field_dependent  = true             # Account for the field, if direct integration is used
relativistic     = true             # Account for relativistic energies
WPI              = true             # Include wave-particle-interactions
adjusted_model   = false            # Use adjusted form of Dentons density model
n                = 1                # Resonance, 1 is the fundamental cyclotron resonance
duct_enhancement = 0.5              # Fraction of density enhancement applied
station          = :sodankyla       # Coordinates of the station simulated for
wave_direction   = -1               # Direction of propagation for the wave, -1 is southward

## Define magnetic field
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

## Define the particle
α_grid = deg2rad.([0.25, 0.5, 0.75, 1.0])
μ_grid = -cos.(α_grid)

# Find the L-shell and equatorial position of the fieldline belonging to the station at the
# requested time
trace_L_shell = find_L_shell(magnetic_field, ref, t)
L = trace_L_shell.L
position_eq = trace_L_shell.position_eq

ϕ_eq = magnetic_field.ϕ
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))
E_grid = range(7e3, 40e3, length=20)

particles = Matrix{ParticleState}(undef, length(E_grid), length(μ_grid))
@showprogress for iα in eachindex(α_grid), (iE, E) in enumerate(E_grid)
    particles[iE, iα] = ParticleState(
        E,
        μ_grid[iα],
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=relativistic,
    )
end

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
lb_low, lb_high = 0.1, 0.5
density = DentonDensity(
    magnetic_field,
    L,
    ϕ_eq;
    duct_enhancement=duct_enhancement,
    adjusted_model=adjusted_model,
)

plasma = PlasmaState(λ_grid, density; lb_low=lb_low, lb_high=lb_high)

## Define the wave
ω_frac_used = 0.3
ωs = [ω_frac_used * minimum(plasma.Ω_e)]

a_Bw = 100e-12
amplitude = ConstantAmplitude(a_Bw)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

Ψ0_grid = collect(range(0.0, 2π, length=5))[1:end-1]
waves = [
    WPI ? WaveState(modulation, amplitude, plasma; Ψ0=Ψ0, wave_direction=wave_direction) : nothing
    for Ψ0 in Ψ0_grid
]


## Calculate TOF and footpoints
TOF = fill(NaN, length(E_grid), length(μ_grid), length(Ψ0_grid))
FP  = fill(NaN, 3, length(E_grid), length(μ_grid), length(Ψ0_grid))

pairs_Eα = CartesianIndices((length(E_grid), length(μ_grid)))
progress = Progress(length(pairs_Eα))
Threads.@threads for k in eachindex(pairs_Eα)
    iE, iα = Tuple(pairs_Eα[k])

    for iΨ in eachindex(Ψ0_grid)
        result = WPI_TOF(
            ωs,
            particles[iE, iα],
            plasma;
            wave=waves[iΨ],
            use_boris=use_boris,
            field_dependent=field_dependent,
            n=n,
            r_start=position_eq,
            get_footpoint=true,
        )
        TOF[iE, iα, iΨ] = result.tof[1]
        FP[:, iE, iα, iΨ] = result.footpoint[:, 1]
    end

    next!(progress)
end
##
params = (
    # method
    use_boris        = use_boris,
    relativistic     = relativistic,
    WPI              = WPI,
    field_dependent  = field_dependent,
    n                = n,
    adjusted_model   = adjusted_model,
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
    ω_frac_used      = ω_frac_used,
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
                       "Particles built at r0 = Dipolar(L, 0, ϕ_eq); r_start = position_eq for every case.",
)

jldsave(name;
    TOF,                          # TOF [E, α, Ψ0]
    FP,                           # FP [xyz, E, α, Ψ0]
    E_grid  = collect(E_grid),
    α_grid,
    μ_grid,
    Ψ0_grid,
    ωs,
    params,
)
## Convert to East-North-Up coordinates with chosen station as reference
reference = LTPReference(station)
fp_enu = fill(NaN, size(FP))
for iΨ in axes(FP, 4), iα in axes(FP, 3), iE in axes(FP, 2)
    fp_enu[:, iE, iα, iΨ] .= FP_to_ENU(FP[:, iE, iα, iΨ], reference, t)
end


##
E_keV = collect(E_grid) ./ 1e3

east = fp_enu[1, :, :, :] ./ 1e3
north = fp_enu[2, :, :, :] ./ 1e3

## Some stats
n_total = length(TOF)
n_no_res = count(isnan, TOF)
n_bounce = count(iszero, TOF)
n_precip = count(x -> isfinite(x) && x > 0, TOF)

## Figure: one axis per pitch angle, stacked
n_rows, n_cols = 2, 2
fig = Figure(size = (1100, 400))
axs = Matrix{Axis}(undef, n_rows, n_cols)

sc = nothing
for iα in eachindex(α_grid)
    r, c = fldmod1(iα, n_cols)
    ax = Axis(fig[r, c];
        xlabel = "East [km]",
        ylabel = "North [km]",
        title  = "α = $(round(rad2deg(α_grid[iα]), digits=2)) deg",
        aspect = DataAspect(),
    )
    axs[r, c] = ax

    # One scatter per phase, coloured by energy
    for iΨ in eachindex(Ψ0_grid)
        sc = scatter!(ax, east[:, iα, iΨ], north[:, iα, iΨ];
            color      = E_keV,
            colormap   = :turbo,
            colorrange = extrema(E_keV),
            markersize = 12,
        )
    end
end

# Same limits in all four panels, then hide the repeated x labels
linkaxes!(vec(axs)...)
for r in 1:n_rows, c in 1:n_cols
    r < n_rows && hidexdecorations!(axs[r, c]; grid=false)
    c > 1 && hideydecorations!(axs[r, c]; grid=false)
end
colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

Colorbar(fig[1:n_rows, n_cols + 1], sc; label = "E [keV]")

config_info = join([
    "Coorinates: $(LTPReference(station).lat_deg)°, $(LTPReference(station).lon_deg)°",
    "Time: $(t)",
    "L = $(round(L, digits=2))",
    "",
    "pdyn = $(pdyn) nPa",
    "Dst = $(dst) nT",
    "By = $(byimf) nT",
    "Bz = $(bzimf) nT",
    "tilt = $(round(rad2deg(magnetic_field.ps), digits=1))°",
    "",
    "ω / Ωₑ = $(ω_frac_used)",
    "Bw = $(round(a_Bw * 1e12, digits=0)) pT",
    "Phases: $(length(Ψ0_grid))",
    "",
    "Adjusted density: $(adjusted_model)",
    "n = $(n)",
    "Duct enhancement: $(duct_enhancement)nₑ₀",
    "",
    "# Particles: $(n_total)",
    "# Precipitated: $(n_precip)",
    "# Bounced: $(n_bounce)",
    "# No resonance: $(n_no_res)" #TODO: Check that this is correct
], "\n")

Label(fig[1:n_rows, n_cols + 2], config_info;
    justification = :left,        # left-aligned lines inside the box
    halign        = :left,
    valign        = :top,
    tellheight    = false,        # the box does not set the height of the rows
    fontsize      = 11,
)

Label(fig[0, 1:n_cols + 2],
    "Particle drift relative to fieldline footpoint of $(station) station";
    fontsize=16, font=:bold, tellwidth=false)

ylims!(axs[1], -0.8, -0.4)


save("src/WPI/test_scripts/footpoint/figures/ductedfootpoints_enu.png", fig)
