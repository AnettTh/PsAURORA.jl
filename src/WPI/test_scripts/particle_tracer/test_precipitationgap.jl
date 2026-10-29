using AURORA
using WGLMakie
using Bonito
using ProgressMeter

## Define a magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
#magnetic_field = DipoleMagneticField()

# Define a particle
μ = -cos(deg2rad(1))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

E_grid = 10e3:1e3:40e3
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

## Define a wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

# Find the resonance latitude
ω_start = wave_frequency(wave, 0.0)
λ_resonance = [resonance_latitude(ω_start, p, wave.plasma)[1] for p in particles]

## Run Boris-mover
results = Vector{Any}(undef, length(E_grid))

progress = Progress(length(E_grid), desc="Energy loop: ")
for (i, p) in enumerate(particles)
    results[i] = boris_mover(p;
    wave=wave,
    n_T=100000,
    resolution=10,
    store_trajectory=true,
    λ_start=λ_resonance[i]
    )
    next!(progress; showvalues=[(:E_keV, round(E_grid[i]/1e3, digits=1))])
end


## Look at footpoints
@load "src/WPI/data/boris_results_L6.5_alpha1.0.jld2" results E_grid λ_resonance L μ ϕ_eq
for (i, r) in enumerate(results)
    @show E_grid[i], r.footpoint
end


## Make figure
fig = Figure(size=(900, 600))
ax  = Axis3(fig[1,1];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    zlabel = "Z [RE]",
    title  = "Ionospheric footpoints",
    aspect = :data
)

# Colormap for energy
colors = cgrad(:viridis, length(E_grid), categorical=true)

# Plot footpoints
for (i, r) in enumerate(results)
    if !isnothing(r.footpoint)
        x_fp, y_fp, z_fp = r.footpoint
        scatter!(ax, [x_fp/RE], [y_fp/RE], [z_fp/RE];
            color      = colors[i],
            markersize = 15,
            marker     = :circle,
            #label      = "E=$(round(E_grid[i]/1e3, digits=1)) keV"
        )
    else
        scatter!(ax, [NaN], [NaN], [NaN];
            color      = :red,
            markersize = 15,
            marker     = :xcross,
            label      = "E=$(round(E_grid[i]/1e3, digits=1)) keV"
        )

    end
end

Colorbar(fig[1,2];
    colormap = :viridis,
    limits   = (minimum(E_grid)/1e3, maximum(E_grid)/1e3),
    label    = "E [keV]"
)

axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
