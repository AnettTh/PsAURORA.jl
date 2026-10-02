using AURORA
using CairoMakie
using Profile
using ProgressMeter


## Define magnetic field
ϕ = 0.0
dipole = DipoleMagneticField()
tsyganenko = TsyganenkoMagneticField(ϕ)

## Define the particles
μ = -cos(deg2rad(3))       # Almost field-aligned
L = 6.0
r0 = Cartesian(L*RE, 0.0, 0.0)

E_grid = range(1e3, 40e3, length=100)

particle_dipole = [ParticleState(E, μ, r0, dipole; relativistic=true) for E in E_grid]
particle_tsyg   = [ParticleState(E, μ, r0, tsyganenko; relativistic=true) for E in E_grid]

## Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
plasma_dipole = PlasmaState(λ_grid, ϕ, DentonDensity(dipole, L, ϕ), dipole, Float64(L))
plasma_tsyg = PlasmaState(λ_grid, ϕ, DentonDensity(tsyganenko, L, ϕ), tsyganenko, Float64(L))

# Define the wave
ωs_dipole = range(plasma_dipole.Ω_e[1]*0.1, plasma_dipole.Ω_e[1]*0.4, length=5)
ωs_tsyg = range(plasma_tsyg.Ω_e[1]*0.1, plasma_tsyg.Ω_e[1]*0.4, length=5)



## Calculate TOF
TOF_dipole = zeros(length(E_grid), length(ωs_dipole))  # [n_E × n_ω]
TOF_tsyg   = zeros(length(E_grid), length(ωs_tsyg  ))


@showprogress for (i, p) in enumerate(particle_dipole)
    TOF_dipole[i, :] = WPI_TOF(ωs_dipole, p, plasma_dipole; field_dependent=true, wave_launch_time=wave_chirp)
end

### Check for bottleneck
#WPI_TOF(ωs_dipole[1], particle_dipole[1], plasma_dipole; field_dependent=true, wave_launch_time=wave_chirp)
#@profview for _ in 1:10
#    WPI_TOF(ωs_dipole[1], particle_dipole[1], plasma_dipole; field_dependent=true, wave_launch_time=wave_chirp)
#end

##
@showprogress for (i, p) in enumerate(particle_tsyg)
    TOF_tsyg[i, :] = WPI_TOF(ωs_tsyg, p, plasma_tsyg; field_dependent=true, wave_launch_time=wave_chirp)
end


## Check for bottleneck
#WPI_TOF(ωs_tsyg[1], particle_tsyg[1], plasma_tsyg; field_dependent=true, wave_launch_time=wave_chirp)
#@profview for _ in 1:10
#    WPI_TOF(ωs_tsyg[1], particle_tsyg[1], plasma_tsyg; field_dependent=true, wave_launch_time=wave_chirp)
#end


## Make Figure
fig = Figure(size=(900,500))
ax = Axis(
    fig[1,1],
    xlabel="time-of-flight [s]",
    ylabel="E [keV]",
    title="Pitch-angle $(round(rad2deg(acos(abs(μ)))))",
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

lines!(ax, [NaN], [NaN]; color=:black, linestyle=:solid, label="Dipole field")
lines!(ax, [NaN], [NaN]; color=:black, linestyle=:dash,  label="Tsyganenko")

Legend(fig[1,2], ax)

##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

#save("src/WPI/test_scripts/energies/tsyg_tof.png", fig)
