using AURORA
using CairoMakie
using Profile
using ProgressMeter
using JLD2

# Running this with E_grid = range(1e3, 40e3, length=50) takes ∼20 minutes
## Define magnetic field
ϕ_eq = deg2rad(240.0)
dipole = DipoleMagneticField()
tsyganenko = TsyganenkoMagneticField(ϕ_eq)

## Define the particles
μ = -cos(deg2rad(2))
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

E_grid = range(1e3, 40e3, length=10)

particle_dipole = Vector{ParticleState}(undef, length(E_grid))
@showprogress for (i, E) in enumerate(E_grid)
    particle_dipole[i] = ParticleState(
        E,
        μ,
        r0,
        dipole,
        ϕ_eq;
        relativistic=true
    )
end
particle_tsyg = Vector{ParticleState}(undef, length(E_grid))
@showprogress for (i, E) in enumerate(E_grid)
    particle_tsyg[i] = ParticleState(
        E,
        μ,
        r0,
        tsyganenko,
        ϕ_eq;
        relativistic=true
    )
end

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=1000)

lb_low, lb_high = 0.1, 0.5
plasma_dipole = PlasmaState(λ_grid, ϕ_eq, DentonDensity(dipole, L, ϕ_eq), dipole, Float64(L); lb_low=lb_low, lb_high=lb_high)
plasma_tsyg = PlasmaState(λ_grid, ϕ_eq, DentonDensity(tsyganenko, L, ϕ_eq), tsyganenko, Float64(L); lb_low=lb_low, lb_high=lb_high)

## Define the wave
ωs_dipole = range(plasma_dipole.ω_lb[1], plasma_dipole.ω_lb[end], length=5)
ωs_tsyg = range(plasma_tsyg.ω_lb[1], plasma_tsyg.ω_lb[end], length=5)

# NOTE: I do not think this is correct anymore
amplitude = ConstantAmplitude(100e-12)                      # 100 pT

modulation_tsyg = Chirp(plasma_tsyg.ω_lb[1], plasma_tsyg.ω_lb[end], 0.2)   # 0.2 s chirp
wave_tsyg = WaveState(modulation_tsyg, amplitude, plasma_tsyg)

modulation_dipole = Chirp(plasma_dipole.ω_lb[1], plasma_dipole.ω_lb[end], 0.2)   # 0.2 s chirp
wave_dipole = WaveState(modulation_dipole, amplitude, plasma_dipole)

## Calculate TOF
TOF_dipole = zeros(length(E_grid), length(ωs_dipole))
TOF_tsyg = zeros(length(E_grid), length(ωs_tsyg))

progress = Progress(length(particle_dipole))
Threads.@threads for i in eachindex(particle_dipole)
    TOF_dipole[i, :] = WPI_TOF(ωs_dipole, particle_dipole[i], plasma_dipole, wave_dipole; use_boris=true, field_dependent=true)
    next!(progress)
end

progress = Progress(length(particle_tsyg))
Threads.@threads for i in eachindex(particle_tsyg)
    TOF_tsyg[i, :] = WPI_TOF(ωs_tsyg, particle_tsyg[i], plasma_tsyg, wave_tsyg; use_boris=true, field_dependent=true)
    next!(progress)
end


## Save the data
#@save "src/WPI/data/TOF_results_L$(L)_alpha$(round(rad2deg(acos(abs(μ))), digits=1)).jld2" TOF_dipole TOF_tsyg E_grid ωs_dipole ωs_tsyg L μ ϕ_eq lb_low lb_high

## Make the figure
fig = Figure(size=(900, 500))
ax = Axis(
    fig[1,1],
    xlabel="TOF [s]",
    ylabel="E [keV]",
    title="Energy as a function of time-of-flight at $(longitude_to_MLT(ϕ_eq; degrees=false)) MLT with pitch-angle $(round(rad2deg(acos(abs(μ))), digits=1))°"
    #yscale=log10
)

colors = [:blue, :red, :green, :orange, :purple]  # one per ω

ω_norm_dipole = range(lb_low, lb_high, length=length(ωs_dipole))
ω_norm_tsyg   = range(lb_low, lb_high, length=length(ωs_tsyg))

for (i, ω) in enumerate(ωs_dipole)
    lines!(ax, TOF_dipole[:, i], E_grid ./ 1e3;
        color     = colors[i],
        linestyle = :solid,
        label     = "ω/Ωₑ = $(round(ω_norm_dipole[i], digits=2))"
    )
end

for (i, ω) in enumerate(ωs_tsyg)
    lines!(ax, TOF_tsyg[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:dash
    )
end

lines!(ax, [NaN], [NaN], label="Dipole field"; color=:black, linestyle=:solid)
lines!(ax, [NaN], [NaN], label="Tsyganenko"; color=:black, linestyle=:dash)

Legend(fig[1,2], ax)

##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

#save("src/WPI/test_scripts/energies/tsyg_tof.png", fig)
save("src/WPI/test_scripts/energies/tsyg_tof_ϕ_$(round(rad2deg(ϕ), digits = 0)).png", fig)
