using AURORA
using Dates
using ProgressMeter
using CairoMakie

## Configuration
duct_enhancement = 0.0
station = :sodankyla
α0 = 0.5
name = "src/WPI/data/11-23-2026/AlphaFootpoint_α0_$(α0).jld2"

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define particle
E = range(7e3, 40e3, length=500)
μ0 = -cos(deg2rad(α0))

particle = [ParticleState(e, μ0, fieldline) for e in E]

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

Ψ0_grid = range(0.0, 2π, length=17)[1:end-1]
waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Frequencies to evaluate
ω_frac = 0.3
ω = ω_frac * minimum(plasma.Ω_e)

nE, nΨ = length(E), length(Ψ0_grid)
α_FP = fill(NaN, nE, nΨ)
FP = fill(NaN, 3, nE, nΨ)
TOF = fill(NaN, nE, nΨ)

# Run every pair of energies and phases
pairs_EΨ = CartesianIndices((nE, nΨ))
progress = Progress(length(pairs_EΨ))
Threads.@threads :greedy for k in eachindex(pairs_EΨ)
    iE, iΨ = Tuple(pairs_EΨ[k])
    # TODO: Run with complete WPI
    result = WPI_TOF(ω, particle[iE], fieldline, plasma; wave=waves[iΨ])
    #result = boris_mover(particle[iE], fieldline, waves[iΨ]; ω=ω)
    if !isnothing(result.footpoint)
        α_FP[iE, iΨ] = result.α_footpoint
        FP[:, iE, iΨ] .= result.footpoint
        TOF[iE, iΨ] = result.tof
    end
    next!(progress)
end

fp_enu = fill(NaN, size(FP))
for iΨ in 1:nΨ, iE in 1:nE
    any(isnan, @view FP[:, iE, iΨ]) && continue
    fp_enu[:, iE, iΨ] .= FP_to_ENU(FP[:, iE, iΨ], ref, t)
end
east_km = fp_enu[1, :, :] ./ 1e3

params = (
    # method
    duct_enhancement = duct_enhancement,
    α0               = α0,
    # station and field
    station          = station,
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
    ω_frac           = ω_frac,
    a_Bw             = amplitude.Bw,                       # [T]
    chirp_duration   = modulation.t_duration,              # [s]
)

jldsave(name;
    TOF,                          # TOF [E, μ, Ψ0]
    FP,                           # FP [xyz, E, μ, Ψ0]
    fp_enu,
    α_FP,
    E_grid=collect(E),
    Ψ0_grid,
    params,
)

E_keV = collect(E) ./ 1e3
E_all = repeat(E_keV, 1, nΨ)
valid = isfinite.(α_FP)
east_lim = maximum(east_km[valid])
TOF_lim = extrema(TOF[valid])

##
fig = Figure(size=(900, 700))
ax1 = Axis(
    fig[1, 1],
    ylabel="α (footpoint) [°]",
    title="Pitch angle at footpoint for precipitating particles at $station $t \nα₀ = $α0"
)

sc1 = scatter!(ax1, E_all[valid], rad2deg.(α_FP[valid]); color=east_km[valid], colormap=:turbo, colorrange=(0, east_lim))
Colorbar(fig[1, 2], sc1; label="East drift [km]")
ylims!(ax1, 0, 90)
xlims!(ax1, 7, 40)

ax2 = Axis(
    fig[2, 1],
    xlabel="E [keV]",
    ylabel="α (footpoint) [°]",
)

sc2 = scatter!(ax2, E_all[valid], rad2deg.(α_FP[valid]); color=TOF[valid], colormap=:viridis, colorrange=TOF_lim)
Colorbar(fig[2, 2], sc2; label="Time of flight [s]")

linkaxes!(ax1, ax2)
hidexdecorations!(ax1; grid=false)
ylims!(ax1, 0, 90)
xlims!(ax1, 7, 40)

rowgap!(fig.layout, 10)

save("src/WPI/test_scripts/footpoint/figures/11232026_AlphaFootpointAndTOF_α$(α0).png", fig)
