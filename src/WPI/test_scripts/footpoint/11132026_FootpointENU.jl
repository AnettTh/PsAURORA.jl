using AURORA
using CairoMakie
using JLD2
using Dates
using ProgressMeter

## Define configuration
name = "src/WPI/data/11-16-2026/FootpointENU_nonducted.jld2"
duct_enhancement = 0.0
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

## Define the particle
α0_grid = [0.25, 0.5, 0.75, 1.0]
μ0_grid = -cos.(deg2rad.(α0_grid))
E_grid = range(10e3, 30e3, length=100)

particles = Matrix{ParticleState}(undef, length(E_grid), length(μ0_grid))
for iμ in eachindex(μ0_grid), (iE, E) in enumerate(E_grid)
    particles[iE, iμ] = ParticleState(E, μ0_grid[iμ], fieldline)
end

## Define the plasma
density = DentonDensity(fieldline; duct_enhancement=duct_enhancement)
plasma = PlasmaState(fieldline, density)

## Make the wave
amplitude = ConstantAmplitude(100e-12)
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

Ψ0_grid = range(0, 3π/2, length=4)
waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Choose frequencies to sample
frac = 0.3
ω = [frac * minimum(plasma.Ω_e)]

## Calculate TOF and footpoints
TOF = fill(NaN, length(E_grid), length(μ0_grid), length(Ψ0_grid))
FP  = fill(NaN, 3, length(E_grid), length(μ0_grid), length(Ψ0_grid))

pairs_Eμ = CartesianIndices((length(E_grid), length(μ0_grid)))
progress = Progress(length(pairs_Eμ))
Threads.@threads for k in eachindex(pairs_Eμ)
    iE, iμ = Tuple(pairs_Eμ[k])

    for iΨ in eachindex(Ψ0_grid)
        result = WPI_TOF(
            ω,
            particles[iE, iμ],
            fieldline,
            plasma;
            wave=waves[iΨ],
            get_footpoint=true,
        )
        TOF[iE, iμ, iΨ] = result.tof[1]
        FP[:, iE, iμ, iΨ] = result.footpoint[:, 1]
    end

    next!(progress)
end
##
params = (
    # method
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
    Ω_eq             = minimum(plasma.Ω_e),                # [rad/s]
    L_Denton_bin     = density.L_Denton_bin,
    R_max            = density.R_max,                      # [m]
    # Wave
    ω_frac           = frac,
    a_Bw             = amplitude.Bw,                       # [T]
    chirp_duration   = modulation.t_duration,              # [s]
)

jldsave(name;
    TOF,                          # TOF [E, μ, Ψ0]
    FP,                           # FP [xyz, E, μ, Ψ0]
    E_grid  = collect(E_grid),
    α0_grid,
    Ψ0_grid,
    params,
)

@load "src/WPI/data/11-16-2026/FootpointENU_nonducted.jld2"

## Convert to East-North-Up coordinates with chosen station as reference
reference = LTPReference(params.station)
fp_enu = fill(NaN, size(FP))
for iΨ in axes(FP, 4), iα in axes(FP, 3), iE in axes(FP, 2)
    fp_enu[:, iE, iα, iΨ] .= FP_to_ENU(FP[:, iE, iα, iΨ], LTPReference(params.station), DateTime(params.time))
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
fig = Figure(size = (900, 400))
axs = Matrix{Axis}(undef, n_rows, n_cols)

scs = nothing
for iα in eachindex(α0_grid)
    r, c = fldmod1(iα, n_cols)
    ax = Axis(fig[r, c];
        xlabel = "East [km]",
        ylabel = "North [km]",
        title  = "α = $(round(α0_grid[iα], digits=2))°",
        aspect = DataAspect(),
    )
    axs[r, c] = ax

    # One scatter per phase, coloured by energy
    for iΨ in eachindex(Ψ0_grid)
        scs = scatter!(ax, east[:, iα, iΨ], north[:, iα, iΨ];
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

Colorbar(fig[1:n_rows, n_cols + 1], scs; label = "E [keV]")

config_info = join([
    "Coorinates: $(LTPReference(params.station).lat_deg)°, $(LTPReference(params.station).lon_deg)°",
    "Time: $(params.time)",
    "",
    "pdyn = $(params.pdyn) nPa",
    "Dst = $(params.dst) nT",
    "By = $(params.byimf) nT",
    "Bz = $(params.bzimf) nT",
    "tilt = $(round(rad2deg(params.ps), digits=1))°",
    "",
    "ω / Ωₑ = $(params.ω_frac)",
    "Bw = $(round(params.a_Bw * 1e12, digits=0)) pT",
    "Phases: $(length(Ψ0_grid))",
    #"",
    #"Duct enhancement: $(duct_enhancement)nₑ₀",
    "",
    "# Particles: $(n_total)",
    "# Precipitated: $(n_precip)",
    "# Bounced: $(n_bounce)",
    "# No resonance: $(n_no_res)"
], "\n")

Label(fig[1:n_rows, n_cols + 2], config_info;
    justification = :left,
    halign        = :left,
    valign        = :top,
    tellheight    = false,
    fontsize      = 11,
)

Label(fig[0, 1:n_cols + 2],
    "Particle drift relative to fieldline footpoint of $(params.station) station";
    fontsize=16, font=:bold, tellwidth=false)

ylims!(axs[1], 4.2, 5.0)


#save("src/WPI/test_scripts/footpoint/figures/11162026_NonDuctedFP_ENU.png", fig)
