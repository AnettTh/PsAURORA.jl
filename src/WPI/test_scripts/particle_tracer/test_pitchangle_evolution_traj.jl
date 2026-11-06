using AURORA
using ProgressMeter
using CairoMakie
using Dates

## Define the magnetic fields
t = DateTime("2020-12-21T03:00:00")
station = :sodankyla
ref = LTPReference(station)

# Storm conditions
byimf = 0.0
bzimf = 2.0
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

## Define a plasma
λ_grid = range(deg2rad(-50.0), deg2rad(50.0), length=200)

density_tsyg = DentonDensity(tsyg, L_tsyg, ϕ_eq; adjusted_model=true, hemisphere=:north)
density_dipole = DentonDensity(dipole, L_dipole, ϕ_eq; adjusted_model=false)

plasma_tsyg = PlasmaState(λ_grid, density_tsyg; lb_low=0.1, lb_high=0.5)
plasma_dipole = PlasmaState(λ_grid, density_dipole; lb_low=0.1, lb_high=0.5)

# Define the waves
modulation_tsyg = Chirp(plasma_tsyg.ω_lb[1], plasma_tsyg.ω_lb[end], 0.2)   # 0.2 s chirp
modulation_dipole = Chirp(plasma_dipole.ω_lb[1], plasma_dipole.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
Ψ0 = 0.0

wave_tsyg = WaveState(modulation_tsyg, amplitude, plasma_tsyg; Ψ0=Ψ0)
wave_dipole = WaveState(modulation_dipole, amplitude, plasma_dipole; Ψ0=Ψ0)


# Find the resonance latitudes
ω_start_tsyg = wave_frequency(wave_tsyg, 0.0)
λ_resonance_tsyg = [resonance_latitude(ω_start_tsyg, p, wave_tsyg.plasma)[1] for p in particles_tsyg]

ω_start_dipole = wave_frequency(wave_dipole, 0.0)
λ_resonance_dipole = [resonance_latitude(ω_start_dipole, p, wave_dipole.plasma)[1] for p in particles_dipole]
## Run Boris-mover (skip particles without resonance)
results_tsyg = Vector{Any}(nothing, length(E_grid))
results_dipole = Vector{Any}(nothing, length(E_grid))

@showprogress for i in eachindex(particles_tsyg)
    if !isnan(λ_resonance_tsyg[i])
        results_tsyg[i] = boris_mover(particles_tsyg[i];
            wave               = wave_tsyg,
            n_T                = 100_000,
            resolution         = 10,
            store_trajectory   = true,
            λ_start            = λ_resonance_tsyg[i],
            stop_when_bouncing = true,
        )
    end
    if !isnan(λ_resonance_dipole[i])
        results_dipole[i] = boris_mover(particles_dipole[i];
            wave               = wave_dipole,
            n_T                = 100_000,
            resolution         = 10,
            store_trajectory   = true,
            λ_start            = λ_resonance_dipole[i],
            stop_when_bouncing = true,
        )
    end
end

## Figure
colors = cgrad(:viridis, length(E_grid), categorical=true)

fig = Figure(size=(900, 700))
ax3 = Axis3(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    zlabel = "Z [RE]",
    aspect = :data,
)


## Earth (drawn first, so the trajectories are layered on top of it)
u = range(0, π,  length=30)
v = range(0, 2π, length=30)
surface!(ax3,
    [sin(ui)*cos(vi) for ui in u, vi in v],
    [sin(ui)*sin(vi) for ui in u, vi in v],
    [cos(ui)         for ui in u, vi in v];
    colormap = :blues, alpha = 0.3, shading = NoShading,
)

## L-bin edges from the Denton density model [RE]
L_edges = [L_dipole]

θ_s = range(0, π,  length=40)
φ_s = range(0, 2π, length=60)

for Lc in L_edges
    xs = [Lc * sin(t) * cos(p) for t in θ_s, p in φ_s]
    ys = [Lc * sin(t) * sin(p) for t in θ_s, p in φ_s]
    zs = [Lc * cos(t)          for t in θ_s, p in φ_s]
    wireframe!(ax3, xs, ys, zs; color = (:gray40, 0.5), linewidth = 0.5)
end

## Trajectories
function plot_trajectories!(ax3, results, field_name, linestyle)
    for (i, r) in enumerate(results)
        isnothing(r) && continue
        n = findlast(>(0), r.time)                 # drop unfilled trailing entries
        isnothing(n) && continue

        lines!(ax3,
            r.position[1:n, 1] ./ RE,
            r.position[1:n, 2] ./ RE,
            r.position[1:n, 3] ./ RE;
            color     = colors[i],
            linewidth = 1.5,
            linestyle = linestyle,
            label     = "E = $(round(E_grid[i] / 1e3, digits=1)) keV, $(field_name)",
        )
    end
end

xlims!(ax3, -8, 1)
ylims!(ax3, -8, 1)
zlims!(ax3, -4, 5)

plot_trajectories!(ax3, results_tsyg,   "Tsyg",   :solid)
plot_trajectories!(ax3, results_dipole, "dipole", :dash)
