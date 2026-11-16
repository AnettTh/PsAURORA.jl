using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Define method
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    t;
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)

## Define the stations fieldline
ref = LTPReference(station)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Make particle
α0 = 2.0
μ0 = -cos(deg2rad(α0))

E_grid = range(7e3, 40e3, length=500)

particles = Vector{ParticleState}(undef, length(E_grid))
for (i, E) in enumerate(E_grid)
    particles[i] = ParticleState(E, μ0, fieldline)
end

## Make the plasma
duct_fractions = 0.0:0.1:0.5
nes = [DentonDensity(fieldline; duct_enhancement=d) for d in duct_fractions]

plasmas = [PlasmaState(fieldline, ne) for ne in nes]

## Make the wave
# NOTE: Is the waves identical? Looks like it...
amplitude = ConstantAmplitude(100e-12)
modulation = Chirp(plasmas[1].ω_lb[1], plasmas[1].ω_lb[end], 0.2)

waves = [WaveState(modulation, amplitude, plasma, fieldline) for plasma in plasmas]

## Choose frequency to simulate for
frac = 0.5
ω = frac * minimum(plasmas[1].Ω_e)

## Calculate TOF
TOF = zeros(length(E_grid), length(plasmas))

progress = Progress(length(particles))
Threads.@threads for i in eachindex(particles)
    for j in eachindex(plasmas)
        TOF[i,j] = WPI_TOF([ω], particles[i], fieldline, plasmas[j]; wave=waves[j]).tof[1]
    end
    next!(progress)
end

println("$(count(x -> isfinite(x) && x > 0, TOF))/$(length(TOF)) precipitaded ($(round(count(x -> isfinite(x) && x > 0, TOF) / length(TOF) * 100, digits=1))%)")

## Make figure
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="TOF [s]",
    ylabel="E [keV]",
    title="Ducted vs non-ducted time-of-flight at $(t)",
)

colors = cgrad(:turbo, length(plasmas), categorical=true)

for j in eachindex(plasmas)
    tof_j = replace(TOF[:, j], 0.0 => NaN)
    scatter!(
        ax,
        tof_j,
        E_grid ./ 1e3;
        color = colors[j],
        label = "nₑ = (1 + $(round(duct_fractions[j], digits=1)))nₑ₀"
    )
end

# Add info to Legend
info = [
    "α = $(round(rad2deg(acos(abs(μ0))), digits=1))°",
    "ω = $(frac)Ωₑ",
]

for i in info
    lines!(ax, [NaN], [NaN]; color=:transparent, label=i)
end

Legend(fig[1,2], ax)

save("src/WPI/test_scripts/ducting/figures/11162026_CompareFrac$(frac)Ωe.png", fig)
