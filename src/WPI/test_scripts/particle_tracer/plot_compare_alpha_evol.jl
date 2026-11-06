using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Define a magnetic field
t = DateTime("2020-12-21T03:00:00")
station = :sodankyla
ref = LTPReference(station)

# Storm conditions
byimf = 0.0
bzimf = -1.6 #2.0
dst = -28.0

tsyg = TsyganenkoMagneticField(ref; time=t, byimf=byimf, bzimf=bzimf, dst=dst)
dipole = DipoleMagneticField()

## Define the particles
μ = -cos(deg2rad(0.5))
L_tsyg = find_L_shell(tsyg, ref, t; r_start=RE+10).L
L_dipole = find_L_shell(dipole, ref, t; r_start=RE+10).L

ϕ_eq = tsyg.ϕ

r0_tsyg = Cartesian(Dipolar(L_tsyg, 0.0, ϕ_eq))
r0_dipole = Cartesian(Dipolar(L_dipole, 0.0, ϕ_eq))

E_grid = 15e3:5e3:40e3 #[1e3; 5e3:5e3:40e3]
particles_tsyg = [
    ParticleState(
        E,
        μ,
        r0_tsyg,
        tsyg,
        ϕ_eq;
        relativistic=true
    )
    for E in E_grid
]

particles_dipole = [
    ParticleState(
        E,
        μ,
        r0_dipole,
        dipole,
        ϕ_eq;
        relativistic=true
    )
    for E in E_grid
]

## Define the plasmas
λ_grid = range(deg2rad(-50.0), deg2rad(50.0), length=200)

density_tsyg = DentonDensity(tsyg, L_tsyg, ϕ_eq; adjusted_model=true, hemisphere=:north)
density_dipole = DentonDensity(dipole, L_dipole, ϕ_eq; adjusted_model=false)

plasma_tsyg = PlasmaState(λ_grid, density_tsyg; lb_low=0.1, lb_high=0.5)
plasma_dipole = PlasmaState(λ_grid, density_dipole; lb_low=0.1, lb_high=0.5)

## Define a wave
modulation_tsyg = Chirp(plasma_tsyg.ω_lb[1], plasma_tsyg.ω_lb[end], 0.2)   # 0.2 s chirp
modulation_dipole = Chirp(plasma_dipole.ω_lb[1], plasma_dipole.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
Ψ0 = 0.0

wave_tsyg = WaveState(modulation_tsyg, amplitude, plasma_tsyg; Ψ0=Ψ0)
wave_dipole = WaveState(modulation_dipole, amplitude, plasma_dipole; Ψ0=Ψ0)

## Find the resonance latitude
ω_start_tsyg = wave_frequency(wave_tsyg, 0.0)
λ_resonance_tsyg = [resonance_latitude(ω_start_tsyg, p, wave_tsyg.plasma)[1] for p in particles_tsyg]

ω_start_dipole = wave_frequency(wave_dipole, 0.0)
λ_resonance_dipole = [resonance_latitude(ω_start_dipole, p, wave_dipole.plasma)[1] for p in particles_dipole]


## Run Boris-mover (skip particles without resonance)
results_tsyg = Vector{Any}(nothing, length(E_grid))
results_dipole = Vector{Any}(nothing, length(E_grid))

@showprogress for i in eachindex(particles_tsyg)
    isnan(λ_resonance_tsyg[i]) && continue
    results_tsyg[i] = boris_mover(particles_tsyg[i];
        wave               = wave_tsyg,
        n_T                = 100_000,
        resolution         = 10,
        store_trajectory   = true,
        λ_start            = λ_resonance_tsyg[i],
        stop_when_bouncing = true,
    )
    isnan(λ_resonance_dipole[i]) && continue
    results_dipole[i] = boris_mover(particles_dipole[i];
        wave               = wave_dipole,
        n_T                = 100_000,
        resolution         = 10,
        store_trajectory   = true,
        λ_start            = λ_resonance_dipole[i],
        stop_when_bouncing = true,
    )
end


## Make colorscale
colors = cgrad(:turbo, length(E_grid), categorical=true)

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

Label(fig[0, :], "Pitch-angle evolution for $(station) at $(t) ($(round(longitude_to_MLT(ϕ_eq; degrees=false), digits=1)) MLT)";
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
    "α₀ = $(round(rad2deg(acos(-μ)), digits=1))°",
    "L_tsyg = $(round(L_tsyg, digits=1))",
    "L_dipole = $(round(L_dipole, digits=1))",
]

for i in info
    lines!(ax_tsyg, [NaN], [NaN]; color=:transparent, label=i)
end

linkaxes!(ax_tsyg, ax_dipole)
hidexdecorations!(ax_dipole; grid=false)
#ylims!(ax_tsyg, -0.5, 15)
#xlims!(ax_tsyg, 0.0, 1.3)
Legend(fig[1:2, 2], ax_tsyg)


#save("src/WPI/test_scripts/particle_tracer/figures/compare_alpha_evolution_bzimf$(bzimf)_zoomed.png", fig)
#save("src/WPI/test_scripts/particle_tracer/figures/compare_alpha_evolution_bzimf$(bzimf).png", fig)
