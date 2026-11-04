using AURORA
using ProgressMeter
using WGLMakie
using Bonito

## Define a magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
#magnetic_field = DipoleMagneticField()

# Define a particle
μ = -cos(deg2rad(2))
L = 6.5
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))

E_grid = [1e3; 10e3; 20e3:5e3:40e3] #[1e3; 5e3:5e3:40e3]
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
density = DentonDensity(magnetic_field, L, ϕ_eq)

plasma = PlasmaState(λ_grid, density; lb_low=0.1, lb_high=0.5)

# Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
Ψ0 = 0.0
wave = WaveState(modulation, amplitude, plasma; Ψ0=Ψ0)


# Find the resonance latitude
ω_start = wave_frequency(wave, 0.0)
λ_resonance = [resonance_latitude(ω_start, p, wave.plasma)[1] for p in particles]

## Run Boris-mover
## Run Boris-mover (skip particles without resonance)
results = Vector{Any}(nothing, length(E_grid))

@showprogress for (i, p) in enumerate(particles)
    isnan(λ_resonance[i]) && continue
    results[i] = boris_mover(p;
        wave             = wave,
        n_T              = 100_000,
        resolution       = 10,
        store_trajectory = true,
        λ_start          = λ_resonance[i],
        stop_when_bouncing = true,
    )
end


## Time window to highlight [s]
t0, t1 = 0.085, 0.175

## Split one result into full path and highlighted window
function split_trajectory(r, t0, t1)
    n = findlast(>(0), r.time)            # last filled entry (time[1] = 0 is the start)
    isnothing(n) && return nothing
    full = 1:n
    in_window = findall(t -> t0 <= t <= t1, r.time[full])
    return (
        full   = full,
        window = in_window,
    )
end

parts = [isnothing(r) ? nothing : split_trajectory(r, t0, t1) for r in results]

## Figure
colors    = cgrad(:viridis, length(E_grid), categorical=true)
thin_full = 1     # plot every 10th point of the full path
thin_win  = 1      # plot every 2nd point of the highlighted window

fig = Figure(size=(1300, 600))
ax3 = Axis3(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    zlabel = "Z [RE]",
    title  = "Trajectories, thick line: t = $(t0) to $(t1) s",
    aspect = :data,
)
ax2 = Axis(fig[1, 2];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "xz-projection",
    aspect = DataAspect(),
)

for (i, p) in enumerate(parts)
    (isnothing(p) || isnothing(results[i])) && continue
    r = results[i]

    # Full trajectory, thin
    idx = p.full[1:thin_full:end]
    xf = r.position[idx, 1] ./ RE
    yf = r.position[idx, 2] ./ RE
    zf = r.position[idx, 3] ./ RE

    lines!(ax3, xf, yf, zf;
        color     = colors[i],
        linewidth = 1,
        label     = "E = $(round(E_grid[i] / 1e3, digits=1)) keV",
    )
    lines!(ax2, xf, zf; color = colors[i], linewidth = 1)

    # Highlighted window, thick, drawn after so it sits on top
    isempty(p.window) && continue
    iw = p.window[1:thin_win:end]
    xw = r.position[iw, 1] ./ RE
    yw = r.position[iw, 2] ./ RE
    zw = r.position[iw, 3] ./ RE

    lines!(ax3, xw, yw, zw; color = colors[i], linewidth = 4)
    lines!(ax2, xw, zw;     color = colors[i], linewidth = 4)
end

## Earth
u = range(0, π,  length=30)
v = range(0, 2π, length=30)
surface!(ax3,
    [sin(ui)*cos(vi) for ui in u, vi in v],
    [sin(ui)*sin(vi) for ui in u, vi in v],
    [cos(ui)         for ui in u, vi in v];
    colormap = :blues, alpha = 0.3, shading = NoShading,
)
θ = range(0, 2π, length=200)
lines!(ax2, cos.(θ), sin.(θ); color = :blue)

## L-bin edges from the Denton density model [RE]
L_edges = [6.6] #[6.6, 7.4, 8.2] #[3.3, 3.7, 4.15, 4.65, 5.2, 5.85, 6.6, 7.4, 8.2]

## Option A: spheres of radius L (what you asked for)
θ_s = range(0, π,  length=40)
φ_s = range(0, 2π, length=60)

for Lc in L_edges
    xs = [Lc * sin(t) * cos(p) for t in θ_s, p in φ_s]
    ys = [Lc * sin(t) * sin(p) for t in θ_s, p in φ_s]
    zs = [Lc * cos(t)          for t in θ_s, p in φ_s]
    wireframe!(ax3, xs, ys, zs; color = (:gray40, 0.5), linewidth = 0.5)
end

## Option B: constant-L surfaces, r = L cos(λ)^2 (the real bin boundaries)
φ_s = range(0, 2π, length=60)

#for Lc in L_edges
#    λ_max = acos(sqrt((RE + z_ionosphere) / (Lc * RE)))
#    λ_s   = range(-λ_max, λ_max, length=40)
#
#    xs = [Lc * cos(l)^3 * cos(p)       for l in λ_s, p in φ_s]
#    ys = [Lc * cos(l)^3 * sin(p)       for l in λ_s, p in φ_s]
#    zs = [Lc * cos(l)^2 * sin(l)       for l in λ_s, p in φ_s]
#    wireframe!(ax3, xs, ys, zs; color = (:orange, 0.5), linewidth = 0.5)
#end

Legend(fig[1, 3], ax3)

## Serve in the browser
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
