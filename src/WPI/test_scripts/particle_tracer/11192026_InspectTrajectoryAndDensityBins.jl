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

tsyg = TsyganenkoMagneticField(t; byimf=byimf, bzimf=bzimf, dst=dst)
dipole = DipoleMagneticField()

## Define station fieldline
fieldline_tsyg = StationFieldLine(tsyg, ref, t)
fieldline_dipole = StationFieldLine(dipole, ref, t)

## Define the particles
E_grid = 15e3:5e3:40e3
μ0 = -cosd(0.5)
particles_tsyg = [ParticleState(E, μ0, fieldline_tsyg) for E in E_grid]
particles_dipole = [ParticleState(E, μ0, fieldline_dipole) for E in E_grid]

## Define a plasma
density_tsyg = DentonDensity(fieldline_tsyg)
density_dipole = DentonDensity(fieldline_dipole)

plasma_tsyg = PlasmaState(fieldline_tsyg, density_tsyg)
plasma_dipole = PlasmaState(fieldline_dipole, density_dipole)

## Define the waves
modulation_tsyg = Chirp(plasma_tsyg.ω_lb[1], plasma_tsyg.ω_lb[end], 0.2)
modulation_dipole = Chirp(plasma_dipole.ω_lb[1], plasma_dipole.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)
Ψ0 = 0.0

wave_tsyg = WaveState(modulation_tsyg, amplitude, plasma_tsyg, fieldline_tsyg; Ψ0=Ψ0)
wave_dipole = WaveState(modulation_dipole, amplitude, plasma_dipole, fieldline_dipole; Ψ0=Ψ0)

## Run Boris-mover (skip particles without resonance)
results_tsyg = Vector{Any}(nothing, length(E_grid))
results_dipole = Vector{Any}(nothing, length(E_grid))

@showprogress for i in eachindex(particles_tsyg)

    results_tsyg[i] = boris_mover(particles_tsyg[i];
        wave               = wave_tsyg,
        n_T                = 100_000,
        resolution         = 10,
        store_trajectory   = true,
        stop_when_bouncing = true,
    )
    results_dipole[i] = boris_mover(particles_dipole[i];
        wave               = wave_dipole,
        n_T                = 100_000,
        resolution         = 10,
        store_trajectory   = true,
        stop_when_bouncing = true,
    )
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
L_edges = [7.4]#[3.3, 3.7, 4.15, 4.65, 5.2, 5.85, 6.6, 7.4, 8.2]

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

##
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
