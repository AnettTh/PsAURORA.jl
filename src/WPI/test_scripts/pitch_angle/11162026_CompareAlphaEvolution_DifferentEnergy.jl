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

E_grid = [7e3; 10e3:10e3:40e3] #[1e3; 5e3:5e3:40e3]
particles_tsyg = [ParticleState(E, μ0, fieldline_tsyg) for E in E_grid]
particles_dipole = [ParticleState(E, μ0, fieldline_dipole) for E in E_grid]

## Define the plasmas
density_tsyg = DentonDensity(fieldline_tsyg)
density_dipole = DentonDensity(fieldline_dipole)

plasma_tsyg = PlasmaState(fieldline_tsyg, density_tsyg)
plasma_dipole = PlasmaState(fieldline_dipole, density_dipole)

## Define a wave
modulation_tsyg = ChirpModulation(plasma_tsyg.ω_lb[1], plasma_tsyg.ω_lb[end], 0.2)
modulation_dipole = ChirpModulation(plasma_dipole.ω_lb[1], plasma_dipole.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)
Ψ0 = 0.0

wave_tsyg = WaveState(modulation_tsyg, amplitude, plasma_tsyg, fieldline_tsyg; Ψ0=Ψ0)
wave_dipole = WaveState(modulation_dipole, amplitude, plasma_dipole, fieldline_dipole; Ψ0=Ψ0)


## Run Boris-mover (skip particles without resonance)
results_tsyg = Vector{Any}(nothing, length(E_grid))
results_dipole = Vector{Any}(nothing, length(E_grid))

@showprogress for i in eachindex(particles_tsyg)
    results_tsyg[i] = boris_mover(
        particles_tsyg[i];
        wave = wave_tsyg,
        store_trajectory = true,
    )
    results_dipole[i] = boris_mover(
        particles_dipole[i];
        wave = wave_dipole,
        store_trajectory = true,
    )
end


## Make colorscale
colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"] #cgrad(:turbo, length(E_grid), categorical=true)

# Plot every nth energy, or specific indices
plot_indices = 1:length(E_grid)

function plot_pitchangle!(ax, results, E_grid, colors)
    for i in eachindex(E_grid)
        isnothing(results[i]) && continue
        r = results[i]
        α = replace(rad2deg.(r.α), 0.0 => NaN)
        lines!(ax, r.time, α;
            color = colors[i],
            label = "$(round(E_grid[i] / 1e3, digits=1)) keV",
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

scatter!(ax_tsyg, [NaN], [NaN], label="Energy", color=:transparent)

plot_pitchangle!(ax_tsyg, results_tsyg, E_grid, colors)
plot_pitchangle!(ax_dipole, results_dipole, E_grid, colors)

## Add info to Legend
info = [
    "bzimf = $(bzimf) nT",
    "dst = $(dst) nT",
    "α₀ = $(round(rad2deg(acos(-μ0)), digits=1))°",
]

for i in info
    lines!(ax_tsyg, [NaN], [NaN]; color=:transparent, label=i)
end

linkaxes!(ax_tsyg, ax_dipole)
hidexdecorations!(ax_dipole; grid=false)
Legend(fig[1:2, 2], ax_tsyg)

#save("src/WPI/test_scripts/pitch_angle/figures/11192026_CompareAlphaEvolution_DifferentEnergy.png", fig)
#
#ylims!(ax_tsyg, -0.1, 15.0)
#xlims!(ax_tsyg, 0.1, 0.4)
#save("src/WPI/test_scripts/pitch_angle/figures/11192026_CompareAlphaEvolution_DifferentEnergy_zoomed.png", fig)
