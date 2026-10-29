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

ne = DentonDensity(magnetic_field, L, ϕ_eq)
ne_ducted = DentonDensity(magnetic_field, L, ϕ_eq; duct_enhancement=0.5)

plasma = PlasmaState(λ_grid, ne)
plasma_ducted = PlasmaState(λ_grid, ne_ducted)

## Make the wave
amp = ConstantAmplitude(100e-12)
mod = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WaveState(mod, amp, plasma)
wave_ducted = WaveState(mod, amp, plasma_ducted)

## Choose frequency to simulate for
frac = 0.4
ω = frac * minimum(plasma.Ω_e)

## Calculate TOF
TOF = zeros(length(E_grid))
TOF_ducted = zeros(length(E_grid))


progress = Progress(length(particles))
Threads.@threads for i in eachindex(particles)
    TOF[i] = WPI_TOF(ω, particles[i], plasma; wave=wave, use_boris=true)
    TOF_ducted[i] = WPI_TOF(ω, particles[i], plasma_ducted; wave=wave_ducted, use_boris=true)
    next!(progress)
end

##
TOF = replace(TOF, 0.0 => NaN)
TOF_ducted = replace(TOF_ducted, 0.0 => NaN)

## Make figure
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="TOF [s]",
    ylabel="E [keV]",
    title="Ducted vs non-ducted time-of-flight",
)

scatter!(ax, TOF, E_grid ./1e3, label="Non-ducted")
scatter!(ax, TOF_ducted, E_grid ./1e3, label="Ducted")

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
