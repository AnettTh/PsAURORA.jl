using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Configuration
duct_enhancement = 0.0
station          = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define particle
E = 30e3
μ0 = -cos(deg2rad(10))
particle = ParticleState(E, μ0, fieldline)

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

Ψ0_grid = range(0, 2π, length=6)[1:end-1]

waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Frequencies to evaluate
ω_frac = 0.3
ω =  ω_frac * minimum(plasma.Ω_e)

## Run Boris-mover
results = Vector{Any}(nothing, length(Ψ0_grid))

@showprogress for i in eachindex(waves)
    # TODO: run with resonance
    #results[i] = boris_mover(particle, fieldline, waves[i]; ω=ω, store_trajectory=true, stop_when_bouncing=false, n_T=10_000)
    results[i] = boris_mover(
        particle;
        wave = waves[i],
        store_trajectory = true,
        stop_when_bouncing = false,
        n_T = 10_000
    )
end

## Latitude and pitch angle of each trajectory
function λ_and_α(r)
    n = findlast(>(0), r.time)
    x, y, z = r.position[1:n, 1], r.position[1:n, 2], r.position[1:n, 3]
    λ = rad2deg.(asin.(z ./ sqrt.(x.^2 .+ y.^2 .+ z.^2)))
    return λ, rad2deg.(r.α[1:n])
end

## Make colorscale
colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

## Figure
fig = Figure(size=(900, 400))
ax  = Axis(fig[1, 1];
    xlabel = "λ (SM) [°]",
    ylabel = "α [°]",
    title  = "Pitch angle evolution along the trajectory for different Ψ₀",
)

scatter!(ax, [NaN], [NaN], label="Ψ₀", color=:transparent)

for (i, r) in enumerate(results)
    λ, α = λ_and_α(r)
    lines!(ax, λ, α; color=colors[i], label="$(round(Ψ0_grid[i]/π, digits=1)) π")
end

# Adiabatic pitch angle along the line, without a wave
α_ad = [begin
            a = pitch_angle(particle.α_eq, particle.B_eq, B)
            isnothing(a) ? NaN : rad2deg(a)
        end for B in fieldline.Bs]

lines!(ax, rad2deg.(fieldline.λs), α_ad; color=:black, linestyle=:dash, label="Loss cone")

# Where the particle started, i.e. the resonance point the mover found
λ_start, _ = λ_and_α(results[1])
vlines!(ax, [λ_start[1]]; linestyle=:dot, color=:black, label="λᵣₑₛ")

axislegend(ax; position=:lt)

ylims!(ax, 0, 90)
xlims!(ax, -70, maximum(rad2deg.(fieldline.λs)))
#save("src/WPI/test_scripts/pitch_angle/figures/11192026_AlphaLatitudeEvolution_HighAlpha.png", fig)
##
#xlims!(ax, 18.5, 19)
#ylims!(ax, 25, 32)
#save("src/WPI/test_scripts/pitch_angle/figures/11192026_AlphaLatitudeEvolution_HighAlpha_zoomed.png", fig)
