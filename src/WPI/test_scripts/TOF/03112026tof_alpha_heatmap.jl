using AURORA
using CairoMakie
using Profile
using ProgressMeter
using JLD2


## Define magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define the particles
α0, α1 = deg2rad(0.01), deg2rad(10)
μ_grid = -cos.(range(α0, α1, length=25))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

E = 30e3

particle = [ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true) for μ in μ_grid]

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
density = DentonDensity(magnetic_field, L, ϕ_eq; duct_enhancement=0.5)
plasma = PlasmaState(λ_grid, density)

# Define the wave
Ω_e_eq = minimum(plasma.Ω_e)
ωs = range(Ω_e_eq*0.1, Ω_e_eq*0.5, length=25)

modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma; Ψ0=0.0)

## Calculate TOF
TOF = zeros(length(μ_grid), length(ωs))
use_boris = true
progress = Progress(length(particle))

Threads.@threads :greedy for i in eachindex(particle)
    TOF[i, :] = WPI_TOF(
        ωs,
        particle[i],
        plasma;
        wave=wave,
        use_boris=use_boris,
        field_dependent=true,
    )
    next!(progress)
end

outfile = "src/WPI/data/11-03-2026/03112026-ducted_TOF_energy_heatmap.jld2"

jldsave(outfile;
    TOF,                              # [α, ω], size (25, 25)
    μ_grid  = collect(μ_grid),        # []
    ωs      = collect(ωs),            # [rad/s]
    ω_frac  = collect(ωs) ./ Ω_e_eq,  # ω/Ωₑ,eq, the axis you plot against
    Ω_e_eq,
    E, L, ϕ_eq,
    Bw      = amplitude.Bw,           # [T]
    Ψ0      = wave.Ψ0,
    t_chirp = modulation.t_duration,
    use_boris,
)

## Make Figure
fig = Figure(size=(700,500))
ax = Axis(
    fig[1,1],
    xlabel="ω/Ωₑ",
    ylabel="α [°]",
    title="Time-of-flight using "
    * if use_boris "boris mover \n" else "integration \n " end
    * "E=$(E / 1e3) keV"
    * ", MLT $(round(longitude_to_MLT(ϕ_eq; degrees=false), digits = 1))"
    * ", L=$(L)",
)

hm = heatmap!(ax, ωs ./ Ω_e_eq, rad2deg.(acos.(-μ_grid)), TOF'; colormap=:turbo, colorrange(2.0, 8.0))
Colorbar(fig[1,2], hm; label="TOF [s]")
##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

save("src/WPI/test_scripts/TOF/figures/03112026_ducted_TOF_alpha_heatmap.png", fig)
