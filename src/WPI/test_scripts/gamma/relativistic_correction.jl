using AURORA
using AURORA; c₀
using CairoMakie


## Make magnetic field
ϕ = deg2rad(120.0)
ϕ_eq = deg2rad(120.0)
magnetic_field = TsyganenkoMagneticField(ϕ)

## Define the particles
μ = -cos(deg2rad(2))       # Almost field-aligned
L = 6.0
r0 = Cartesian(L*RE, 0.0, 0.0)

E_grid = range(1e3, 40e3, length=100)

p_relativistic = [ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true) for E in E_grid]
p_classical = [ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=false) for E in E_grid]

## Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
plasma = PlasmaState(λ_grid, ϕ, DentonDensity(magnetic_field, L, ϕ), magnetic_field, Float64(L))


## Define the wave
ωs = range(plasma.Ω_e[1]*0.1, plasma.Ω_e[1]*0.4, length=5)

## Calculate TOF
TOF_matrix_relativistic = zeros(length(E_grid), length(ωs))  # [n_E × n_ω]
TOF_simple_relativistic = zeros(length(E_grid), length(ωs))

for (i, p) in enumerate(p_relativistic)
    TOF_matrix_relativistic[i, :] = WPI_TOF(ωs, p, plasma; field_dependent=true, wave_launch_time=wave_chirp)
    TOF_simple_relativistic[i, :]  = WPI_TOF(ωs, p, plasma; field_dependent=false, wave_launch_time=wave_chirp)
end

TOF_matrix_classical = zeros(length(E_grid), length(ωs))  # [n_E × n_ω]
TOF_simple_classical = zeros(length(E_grid), length(ωs))

for (i, p) in enumerate(p_classical)
    TOF_matrix_classical[i, :] = WPI_TOF(ωs, p, plasma; field_dependent=true, wave_launch_time=wave_chirp)
    TOF_simple_classical[i, :]  = WPI_TOF(ωs, p, plasma; field_dependent=false, wave_launch_time=wave_chirp)
end


##
fig = Figure(size=(1400, 500))
Label(fig[0, 1:2],
    "Energy as a function of `time-of-flight`,  " *
    "α=$(round(rad2deg(acos(abs(μ))), digits=1))°, ϕ=$(round(rad2deg(ϕ), digits=1))°";
    fontsize=16,
    font=:bold,
    tellwidth=false,
    halign=:center)

# Panel 1: Simple mode (field-independent), both classical and relativistic
ax1 = Axis(fig[1, 1],
    xlabel = "time-of-flight [s]",
    ylabel = "E [keV]",
    title  = "Field-independent"
)

# Panel 2: Field-dependent mode, both classical and relativistic
ax2 = Axis(fig[1, 2],
    xlabel = "time-of-flight [s]",
    ylabel = "E [keV]",
    title  = "Field-dependent"
)

colors = [:blue, :red, :green, :orange, :purple]  # one per ω

for (i, ω) in enumerate(ωs)
    # Panel 1: simple mode
    lines!(ax1, TOF_simple_classical[:, i],    E_grid ./ 1e3;
        color=colors[i], linestyle=:solid)
    lines!(ax1, TOF_simple_relativistic[:, i], E_grid ./ 1e3;
        color=colors[i], linestyle=:dash)

    # Panel 2: field-dependent mode
    lines!(ax2, TOF_matrix_classical[:, i],    E_grid ./ 1e3;
        color=colors[i], linestyle=:solid,
        label="ω = $(round((ω * 2π)/1e3, digits=1)) kHz")
    lines!(ax2, TOF_matrix_relativistic[:, i], E_grid ./ 1e3;
        color=colors[i], linestyle=:dash)
end

# Dummy lines for legend
lines!(ax2, [NaN], [NaN]; color=:black, linestyle=:solid, label="Classical")
lines!(ax2, [NaN], [NaN]; color=:black, linestyle=:dash,  label="Relativistic")

Legend(fig[1, 3], ax2)

##
#save("src/WPI/test_scripts/gamma/rel_vs_classic.png", fig)
#save("src/WPI/test_scripts/gamma/tsyg_rel_vs_classic.png", fig)
#save("src/WPI/test_scripts/gamma/tsyg_rel_vs_classic_MLT4.png", fig)
