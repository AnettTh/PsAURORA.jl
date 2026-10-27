using AURORA
using WGLMakie
using Bonito
using ProgressMeter

## Define a magnetic field
ϕ_eq = deg2rad(240.0)
#magnetic_field = TsyganenkoMagneticField(ϕ_eq)
magnetic_field = DipoleMagneticField()

# Define a particle
μ = -cos(deg2rad(1))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

E_grid = 30e3:5e3:40e3
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
ne_model = DentonDensity(magnetic_field, L, ϕ_eq)

plasma = PlasmaState(λ_grid, ϕ_eq, ne_model, magnetic_field, L; lb_low=0.1, lb_high=0.5)

# Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

# Find the resonance latitude
ω_start = wave_frequency(wave, 0.0)
λ_resonance = [resonance_latitude(ω_start, p, wave.plasma)[1] for p in particles]

## Run Boris-mover
results = Vector{Any}(undef, length(E_grid))

@showprogress for (i, p) in enumerate(particles)
    results[i] = boris_mover(p;
    wave=wave,
    n_T=100000,
    resolution=10,
    store_trajectory=true,
    λ_start=λ_resonance[i]
    )
end


##
## Trajectory plot
fig = Figure(size=(800, 600))
ax  = Axis(fig[1,1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Particle trajectories (xz-plane)",
    aspect = DataAspect()
)

# Earth
θ_earth = range(0, 2π, length=200)
lines!(ax, cos.(θ_earth), sin.(θ_earth); color=:blue, linewidth=2, label="Earth")

# Field line
r0_field  = Cartesian(Spherical(particles[1].L, 0.0, particles[1].ϕ_eq))
_, trace  = find_R_max(magnetic_field, r0_field; store_trace=true)
lines!(ax, trace.xs./RE, trace.zs./RE; color=:gray, linestyle=:dash, label="Field line")

# Resonance points
colors = cgrad(:viridis, length(E_grid), categorical=true)
for (i, (r, λ)) in enumerate(zip(results, λ_resonance))
    if !isnan(λ)
        r_res = Cartesian(Spherical(particles[i].L, λ, particles[i].ϕ_eq))
        scatter!(ax, [r_res.x/RE], [r_res.z/RE];
            color=colors[i], markersize=15, marker=:diamond)
    end
end

# Trajectories
for (i, r) in enumerate(results)
    lines!(ax, r.position[:, 1]./RE, r.position[:, 3]./RE;
        color = colors[i],
        label = "E=$(round(E_grid[i]/1e3, digits=1)) keV"
    )
end

axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)

## Pitch-angle plot
fig = Figure(size=(800, 400))
ax  = Axis(fig[1,1];
    xlabel = "t [s]",
    ylabel = "α [°]",
    title  = "Pitch angle evolution"
)

colors = cgrad(:viridis, length(E_grid), categorical=true)

# Plot every nth energy, or specific indices
plot_indices = 1:5:length(E_grid)   # every 5th energy

for i in plot_indices
    r = results[i]
    lines!(ax, r.time, rad2deg.(r.α);
        color = colors[i],
        label = "E=$(round(E_grid[i]/1e3, digits=1)) keV"
    )
end

# Add loss cone line
hlines!(ax, [rad2deg(particles[1].α_lc)];
    color=:red, linestyle=:dash, label="Loss cone")

axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)

##
fig = Figure(size=(800, 400))
ax  = Axis(fig[1,1];
    xlabel = "t [s]",
    ylabel = "α [°]",
    title  = "Pitch angle evolution"
)

lines!(ax, result.time, rad2deg.(result.α))

# Add loss cone line
hlines!(ax, [rad2deg(particle.α_lc)]; color=:red, linestyle=:dash, label="Loss cone")

axislegend(ax)
xlims!(ax, 0, result.time[end])

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
