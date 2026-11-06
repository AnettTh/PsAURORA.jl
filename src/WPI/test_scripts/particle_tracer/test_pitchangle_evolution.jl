using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Define a magnetic field
t = DateTime("2020-12-21T03:00:00")
ref = LTPReference(:sodankyla)

# Storm conditions
byimf = 0.0
bzimf = 2.0
dst = -28.0

#magnetic_field = TsyganenkoMagneticField(ref; time=t, byimf=byimf, bzimf=bzimf, dst=dst)
#ϕ_eq = magnetic_field.ϕ

magnetic_field = DipoleMagneticField()
ϕ_eq = deg2rad(269.6)

## Define the particles
μ = -cos(deg2rad(0.5))
L = find_L_shell(magnetic_field, ref, t; r_start=RE+10).L
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))

E_grid = 15e3:5e3:40e3 #[1e3; 5e3:5e3:40e3]
particles = [
    ParticleState(
        E,
        μ,
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=true
    )
    for E in E_grid
]

## Define a plasma
λ_grid = range(deg2rad(-50.0), deg2rad(50.0), length=200)
density = DentonDensity(magnetic_field, L, ϕ_eq; adjusted_model=true, hemisphere=:north)

plasma = PlasmaState(λ_grid, density; lb_low=0.1, lb_high=0.5)

## Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
Ψ0 = 0.0
wave = WaveState(modulation, amplitude, plasma; Ψ0=Ψ0)


## Find the resonance latitude
ω_start = wave_frequency(wave, 0.0)
λ_resonance = [resonance_latitude(ω_start, p, wave.plasma)[1] for p in particles]

## Run Boris-mover (skip particles without resonance)
results = Vector{Any}(nothing, length(E_grid))

@showprogress for (i, p) in enumerate(particles)
    isnan(λ_resonance[i]) && continue
    results[i] = boris_mover(p;
        wave               = wave,
        n_T                = 100_000,
        resolution         = 10,
        store_trajectory   = true,
        λ_start            = λ_resonance[i],
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
            label = "E=$(round(E_grid[i] / 1e3, digits=1)) keV",
        )
    end
end

## Pitch-angle plot
fig = Figure(size=(900, 600)) #1200, 450))

Label(fig[0, :], "Pitch-angle evolution";
    fontsize  = 18,
    font      = :bold,
    tellwidth = false,
)

ax_full = Axis(fig[1, 1], xlabel = "t [s]", ylabel = "α [°]")
ax_zoom = Axis(fig[2, 1], xlabel = "t [s]", ylabel = "α [°]", title = "Zoom")

plot_pitchangle!(ax_full, results, E_grid, colors)
plot_pitchangle!(ax_zoom, results, E_grid, colors)


xlims!(ax_zoom, 0.14, 0.50)
ylims!(ax_zoom, -0.5, 10.0)

# Add info to Legend
info = [
    "α₀ = $(round(rad2deg(acos(-μ)), digits=1))°",
    "L = $(round(L, digits=1))",
    "$(round(longitude_to_MLT(ϕ_eq; degrees=false), digits=1)) MLT",
]

for i in info
    lines!(ax_full, [NaN], [NaN]; color=:transparent, label=i)
end


Legend(fig[1:2, 2], ax_full)
