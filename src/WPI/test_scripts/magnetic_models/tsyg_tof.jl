using AURORA
using CairoMakie
using Profile
using ProgressMeter
using JLD2


## Define magnetic field
ϕ_eq = deg2rad(240.0)
dipole = DipoleMagneticField()
tsyganenko = TsyganenkoMagneticField(ϕ_eq)

## Define the particles
μ = -cos(deg2rad(2))       # Almost field-aligned
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

##
E_grid = range(10e3, 30e3, length=20)

particle_dipole = [
    ParticleState(
        E,
        μ,
        r0,
        dipole,
        ϕ_eq;
        relativistic=true,
        precomputed_magnetic_field=false
    )
    for E in E_grid
]
particle_tsyg = [
    ParticleState(
        E,
        μ,
        r0,
        tsyganenko,
        ϕ_eq;
        relativistic=true,
        precomputed_magnetic_field=false
    )
    for E in E_grid
]

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=100)#, length=500)
plasma_dipole = PlasmaState(λ_grid, ϕ_eq, DentonDensity(dipole, L, ϕ_eq), dipole, Float64(L))
plasma_tsyg = PlasmaState(λ_grid, ϕ_eq, DentonDensity(tsyganenko, L, ϕ_eq), tsyganenko, Float64(L))

# Define the wave
Ω_e_eq_dipole = minimum(plasma_dipole.Ω_e)
Ω_e_eq_tsyg = minimum(plasma_tsyg.Ω_e)

## TODO: Figure out why this won't run for 0.49*Ω_e
ωs_dipole = range(Ω_e_eq_dipole*0.1, Ω_e_eq_dipole*0.3, length=5)
ωs_tsyg = range(Ω_e_eq_tsyg*0.1, Ω_e_eq_tsyg*0.3, length=5)


## Calculate TOF
TOF_dipole = zeros(length(E_grid), length(ωs_dipole))  # [n_E × n_ω]
TOF_tsyg   = zeros(length(E_grid), length(ωs_tsyg  ))


##
progress = Progress(length(particle_dipole))
Threads.@threads for i in eachindex(particle_dipole)
    TOF_dipole[i, :] = WPI_TOF(ωs_dipole, particle_dipole[i], plasma_dipole; use_boris=true, field_dependent=true, wave_launch_time=wave_chirp)
    next!(progress)
end
##
progress = Progress(length(particle_tsyg))
Threads.@threads for i in eachindex(particle_tsyg)
    TOF_tsyg[i, :] = WPI_TOF(ωs_tsyg, particle_tsyg[i], plasma_tsyg; use_boris=true, field_dependent=true, wave_launch_time=wave_chirp)
    next!(progress)
end

##
#@save "src/WPI/data/TOF_results_L$(L)_alpha$(round(rad2deg(acos(abs(μ))), digits=1)).jld2" begin
#    # Parameters
#    L         = L
#    μ         = μ
#    ϕ_eq      = ϕ_eq
#    E_grid    = collect(E_grid)
#    ωs_dipole = collect(ωs_dipole)
#    ωs_tsyg   = collect(ωs_tsyg)
#
#    # Results
#    TOF_dipole = TOF_dipole
#    TOF_tsyg   = TOF_tsyg
#
#    # Plasma parameters for reference
#    Ω_e_eq_dipole = Ω_e_eq_dipole
#    Ω_e_eq_tsyg   = Ω_e_eq_tsyg
#end

## Make Figure
fig = Figure(size=(900,500))
ax = Axis(
    fig[1,1],
    xlabel="time-of-flight [s]",
    ylabel="E [keV]",
    title="Energy as a function of `time-of-flight`, "
    * "α = $(round(rad2deg(acos(abs(μ))), digits=1)), ϕ = $(round(rad2deg(ϕ_eq), digits = 1))",
    #yscale=log10
)

colors = [:blue, :red, :green, :orange, :purple]  # one per ω

for (i, ω) in enumerate(ωs_dipole)
    lines!(ax, TOF_dipole[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:solid,
        label="ω = $(round(ω * 2π/1e3, digits=1)) kHz"
    )
end
for (i, ω) in enumerate(ωs_tsyg)
    lines!(ax, TOF_tsyg[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:dash,
        label="ω = $(round(ω * 2π/1e3, digits=1)) kHz"
    )
end

lines!(ax, [NaN], [NaN], label="Dipole field"; color=:black, linestyle=:solid)
lines!(ax, [NaN], [NaN], label="Tsyganenko"; color=:black, linestyle=:dash)

Legend(fig[1,2], ax)

##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

#save("src/WPI/test_scripts/energies/tsyg_tof.png", fig)
#save("src/WPI/test_scripts/energies/tsyg_tof_ϕ_$(round(rad2deg(ϕ), digits = 0)).png", fig)
