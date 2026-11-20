using AURORA
using WGLMakie
using Bonito
using ProgressMeter


##
t = DateTime("2020-12-21T02:00:00")
station = :sodankyla

## Make magnetic field
magnetic_field = TsyganenkoMagneticField(t)
#magnetic_field = DipoleMagneticField()
fieldline = StationFieldLine(magnetic_field, LTPReference(station), t)

## Define a particle
μ0 = -cosd(0.5)
E_grid = 30e3:5e3:40e3
particles = [ParticleState(E, μ0, fieldline) for E in E_grid]

## Define a plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define a wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma, fieldline)

# Find the resonance latitude
ω_start = wave_frequency(wave, 0.0)
resonances = [point_of_resonance(ω_start, p, fieldline, plasma, wave) for p in particles]

## Run Boris-mover, only for particles that resonate
results = Vector{Any}(nothing, length(E_grid))

@showprogress for (i, p) in enumerate(particles)
    isnothing(resonances[i]) && continue
    results[i] = boris_mover(p;
        wave=wave,
        r_resonance=resonances[i].r,
        n_T=100_000,
        resolution=10,
        store_trajectory=true,
    )
end

# Drop the unfilled trailing entries of a stored trajectory
function filled(r)
    n = findlast(>(0), r.time)
    return (position=r.position[1:n, :], time=r.time[1:n], α=r.α[1:n])
end

# Distance along the meridian through the magnetic equator, so the line and the
# trajectories can be drawn in one plane
ϕ_eq = fieldline.ϕ_eq
meridian(x, y) = @. x * cos(ϕ_eq) + y * sin(ϕ_eq)

colors = cgrad(:viridis, length(E_grid), categorical=true)

## Trajectory plot
fig = Figure(size=(800, 600))
ax  = Axis(fig[1,1];
    xlabel = "Distance along the equator meridian [RE]",
    ylabel = "Z (SM) [RE]",
    title  = "Particle trajectories (meridional projection)",
    aspect = DataAspect()
)

# Earth
θ_earth = range(0, 2π, length=200)
lines!(ax, cos.(θ_earth), sin.(θ_earth); color=:blue, linewidth=2, label="Earth")

# Field line, straight from the struct
lines!(ax, meridian(fieldline.xs, fieldline.ys) ./ RE, fieldline.zs ./ RE;
    color=:gray, linestyle=:dash, label="Field line")

# Resonance points
for (i, res) in enumerate(resonances)
    isnothing(res) && continue
    scatter!(ax, [meridian(res.r.x, res.r.y) / RE], [res.r.z / RE];
        color=colors[i], markersize=15, marker=:diamond)
end

# Trajectories
for (i, r) in enumerate(results)
    isnothing(r) && continue
    f = filled(r)
    lines!(ax, meridian(f.position[:, 1], f.position[:, 2]) ./ RE, f.position[:, 3] ./ RE;
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

for (i, r) in enumerate(results)
    isnothing(r) && continue
    f = filled(r)
    lines!(ax, f.time, rad2deg.(f.α);
        color = colors[i],
        label = "E=$(round(E_grid[i]/1e3, digits=1)) keV"
    )
end

# Loss cone for the hemisphere the particles move towards (μ0 < 0 is northward)
hlines!(ax, [losscone_angle(fieldline; hemisphere=:north, degrees=true)];
    color=:red, linestyle=:dash, label="Loss cone")

axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)

## Single particle: pitch angle over its own time range
i_one = findfirst(!isnothing, results)
isnothing(i_one) && error("No particle resonated, so there is nothing to plot")
f = filled(results[i_one])

fig = Figure(size=(800, 400))
ax  = Axis(fig[1,1];
    xlabel = "t [s]",
    ylabel = "α [°]",
    title  = "Pitch angle evolution, E = $(round(E_grid[i_one]/1e3, digits=1)) keV"
)

lines!(ax, f.time, rad2deg.(f.α))
hlines!(ax, [losscone_angle(fieldline; hemisphere=:north, degrees=true)];
    color=:red, linestyle=:dash, label="Loss cone")

axislegend(ax)
xlims!(ax, 0, f.time[end])

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
