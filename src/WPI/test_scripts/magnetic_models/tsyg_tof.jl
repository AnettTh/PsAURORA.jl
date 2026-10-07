using AURORA
using CairoMakie
using Profile
using ProgressMeter
using JLD2


## Define magnetic field
ϕ = deg2rad(240.0)
ϕ_eq = deg2rad(240.0)
dipole = DipoleMagneticField()
tsyganenko = TsyganenkoMagneticField(ϕ)

## Define the particles
μ = -cos(deg2rad(2))       # Almost field-aligned
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

E_grid = range(10e3, 30e3, length=50)

##
particle_dipole = [
    ParticleState(
        E,
        μ,
        r0,
        dipole,
        ϕ_eq;
        relativistic=true,
        precomputed_magnetic_field=true
    )
    for E in E_grid
]
particle_tsyg = [ParticleState(E, μ, r0, tsyganenko, ϕ_eq; relativistic=true, precomputed_magnetic_field=true) for E in E_grid]

## Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
plasma_dipole = PlasmaState(λ_grid, ϕ_eq, DentonDensity(dipole, L, ϕ_eq), dipole, Float64(L))
plasma_tsyg = PlasmaState(λ_grid, ϕ_eq, DentonDensity(tsyganenko, L, ϕ_eq), tsyganenko, Float64(L))

# Define the wave
ωs_dipole = range(plasma_dipole.Ω_e[1]*0.1, plasma_dipole.Ω_e[1]*0.4, length=5)
ωs_tsyg = range(plasma_tsyg.Ω_e[1]*0.1, plasma_tsyg.Ω_e[1]*0.4, length=5)



## Calculate TOF
TOF_dipole = zeros(length(E_grid), length(ωs_dipole))  # [n_E × n_ω]
TOF_tsyg   = zeros(length(E_grid), length(ωs_tsyg  ))



@showprogress for (i, p) in enumerate(particle_dipole)
    TOF_dipole[i, :] = WPI_TOF(ωs_dipole, p, plasma_dipole; use_boris=true, field_dependent=true, wave_launch_time=wave_chirp)
end

## Check for bottleneck
#WPI_TOF(ωs_dipole[1], particle_dipole[1], plasma_dipole; use_boris=true, field_dependent=true, wave_launch_time=wave_chirp)
#@profview for _ in 1:10
#    WPI_TOF(ωs_dipole[1], particle_dipole[1], plasma_dipole; use_boris=true, field_dependent=true, wave_launch_time=wave_chirp)
#end

##
@showprogress for (i, p) in enumerate(particle_tsyg)
    TOF_tsyg[i, :] = WPI_TOF(ωs_tsyg, p, plasma_tsyg; use_boris=true, field_dependent=true, wave_launch_time=wave_chirp)
end


## Check for bottleneck
#WPI_TOF(ωs_tsyg[1], particle_tsyg[1], plasma_tsyg; field_dependent=true, wave_launch_time=wave_chirp)
#@profview for _ in 1:10
#    WPI_TOF(ωs_tsyg[1], particle_tsyg[1], plasma_tsyg; field_dependent=true, wave_launch_time=wave_chirp)
#end

##
@save "src/WPI/data/TOF_results.jld2" TOF_dipole TOF_tsyg E_grid ωs_dipole ωs_tsyg


## Make Figure
fig = Figure(size=(900,500))
ax = Axis(
    fig[1,1],
    xlabel="time-of-flight [s]",
    ylabel="E [keV]",
    title="Energy as a function of `time-of-flight`, "
    * "α = $(round(rad2deg(acos(abs(μ))), digits=1)), ϕ = $(round(rad2deg(ϕ), digits = 1))",
    #yscale=log10
)

#ax.yticks = [1, 10, 100, 1000]
#ax.ytickformat = values -> ["$(Int(v))" for v in values]

#ω_norm_dipole = (ωs_dipole .- minimum(ωs_dipole)) ./ (maximum(ωs_dipole) - minimum(ωs_dipole))
#ω_norm_tsyg   = (ωs_tsyg   .- minimum(ωs_tsyg  )) ./ (maximum(ωs_tsyg  ) - minimum(ωs_tsyg  ))


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
