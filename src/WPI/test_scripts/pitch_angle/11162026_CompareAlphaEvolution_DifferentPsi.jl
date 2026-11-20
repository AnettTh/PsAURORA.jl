using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Define configuration
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

tsyg = TsyganenkoMagneticField(
    t;
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)
dipole = DipoleMagneticField()

## Define the stations fieldline
ref = LTPReference(station)
fieldline_tsyg = StationFieldLine(tsyg, ref, t)
fieldline_dipole = StationFieldLine(dipole, ref, t)

## Define the particles
α0 = 0.5
μ0 = -cos(deg2rad(α0))

E = 15e3
particle_tsyg = ParticleState(E, μ0, fieldline_tsyg)
particle_dipole = ParticleState(E, μ0, fieldline_dipole)

## Define the plasmas
density_tsyg = DentonDensity(fieldline_tsyg)
density_dipole = DentonDensity(fieldline_dipole)

plasma_tsyg = PlasmaState(fieldline_tsyg, density_tsyg)
plasma_dipole = PlasmaState(fieldline_dipole, density_dipole)

## Define a wave
modulation_tsyg = ChirpModulation(plasma_tsyg.ω_lb[1], plasma_tsyg.ω_lb[end], 0.2)
modulation_dipole = ChirpModulation(plasma_dipole.ω_lb[1], plasma_dipole.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

Ψ0_grid = range(0, 2π, length=6)[1:end-1]

waves_tsyg = [WaveState(modulation_tsyg, amplitude, plasma_tsyg, fieldline_tsyg; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]
waves_dipole = [WaveState(modulation_dipole, amplitude, plasma_dipole, fieldline_dipole; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Run Boris-mover (skip particles without resonance)
results_tsyg = Vector{Any}(nothing, length(Ψ0_grid))
results_dipole = Vector{Any}(nothing, length(Ψ0_grid))

@showprogress for i in eachindex(waves_tsyg)
    results_tsyg[i] = boris_mover(
        particle_tsyg;
        wave = waves_tsyg[i],
        store_trajectory = true,
    )
    results_dipole[i] = boris_mover(
        particle_dipole;
        wave = waves_dipole[i],
        store_trajectory = true,
    )
end


## Make colorscale
colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

# Plot every nth energy, or specific indices
plot_indices = 1:length(Ψ0_grid)

function plot_pitchangle!(ax, results, Ψ0_grid, colors)
    for i in eachindex(Ψ0_grid)
        isnothing(results[i]) && continue
        r = results[i]
        α = replace(rad2deg.(r.α), 0.0 => NaN)
        lines!(ax, r.time, α;
            color = colors[i],
            label = "$(round(Ψ0_grid[i]/π, digits=1))π",
        )
    end
end

## Pitch-angle plot
fig = Figure(size=(900, 600)) #1200, 450))

Label(fig[0, :], "Pitch-angle evolution for $(station) at $(t)";
    fontsize  = 18,
    font      = :bold,
    tellwidth = false,
)

ax_tsyg = Axis(fig[2, 1], xlabel = "t [s]", ylabel = "α [°]", title = "Tsyganenko")
ax_dipole = Axis(fig[1, 1], xlabel = "t [s]", ylabel = "α [°]", title = "Dipole")

scatter!(ax_tsyg, [NaN], [NaN], label="Ψ₀", color=:transparent)

plot_pitchangle!(ax_tsyg, results_tsyg, Ψ0_grid, colors)
plot_pitchangle!(ax_dipole, results_dipole, Ψ0_grid, colors)

## Add info to Legend
info = [
    "bzimf = $(bzimf) nT",
    "dst = $(dst) nT",
    "α₀ = $(round(rad2deg(acos(-μ0)), digits=1))°",
    "E = $(E / 1e3) keV"
]

for i in info
    lines!(ax_tsyg, [NaN], [NaN]; color=:transparent, label=i)
end

linkaxes!(ax_tsyg, ax_dipole)
hidexdecorations!(ax_dipole; grid=false)
Legend(fig[1:2, 2], ax_tsyg)


#save("src/WPI/test_scripts/pitch_angle/figures/11192026_CompareAlphaEvolution_DifferentPsi.png", fig)
#
#ylims!(ax_tsyg, -0.1, 15.0)
#xlims!(ax_tsyg, 0.2, 0.3)
#save("src/WPI/test_scripts/pitch_angle/figures/11192026_CompareAlphaEvolution_DifferentPsi_zoomed.png", fig)
