using AURORA
using CairoMakie
using ProgressMeter

## Make the magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Make particle
μ = -cos(deg2rad(2))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

E_grid = range(10e3, 40e3, length=50)

particles = Vector{ParticleState}(undef, length(E_grid))
@showprogress for (i, E) in enumerate(E_grid)
    particles[i] = ParticleState(
        E,
        μ,
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=true,
    )
end

## Make the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
duct_fractions = 0.0:0.1:0.5

nes = [
    DentonDensity(magnetic_field, L, ϕ_eq; duct_enhancement=d)
    for d in duct_fractions
]

plasmas = [PlasmaState(λ_grid, ne) for ne in nes]

## Make the wave
amp = ConstantAmplitude(100e-12)
mod = Chirp(plasmas[1].ω_lb[1], plasmas[1].ω_lb[end], 0.2)

wave = [WaveState(mod, amp, plasma) for plasma in plasmas]

## Choose frequency to simulate for
frac = 0.1
ω = frac * minimum(plasmas[1].Ω_e)

## Calculate TOF
TOF = zeros(length(E_grid), length(plasmas))

progress = Progress(length(particles))
Threads.@threads for i in eachindex(particles)
    for j in eachindex(plasmas)
        TOF[i,j] = WPI_TOF([ω], particles[i], plasmas[j]; wave=wave[j], use_boris=true)[1]
    end
    next!(progress)
end


## Make figure
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="TOF [s]",
    ylabel="E [keV]",
    title="Ducted vs non-ducted time-of-flight",
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
    "α = $(round(rad2deg(acos(abs(μ))), digits=1))°",
    "L = $L",
    "ω = $(frac)Ωₑ",
    "$(round(longitude_to_MLT(ϕ_eq; degrees=false), digits=1)) MLT",
]

for i in info
    lines!(ax, [NaN], [NaN]; color=:transparent, label=i)
end

Legend(fig[1,2], ax)
