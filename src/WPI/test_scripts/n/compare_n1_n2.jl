using AURORA
using AURORA; c₀
using CairoMakie

## Define the particles
μ = -cos(deg2rad(3))       # Almost field-aligned
L = 6.5
r0 = [L*RE, 0.0, 0.0]

E_grid = range(1e3, 40e3, length=100)

particles = [ParticleState(E, μ, r0, dipole_field; relativistic=true) for E in E_grid]

## Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
n_e0 = 1.8e7        # 18/cc, from Hsieh 2022

plasma = PlasmaState(λ_grid, n_e0, dipole_field, L)


## Define the wave
Ωe_frac = range(0.1, 0.5, length=5)
ωs = plasma.Ω_e[1] .* Ωe_frac

##
TOF_n1 = zeros(length(E_grid), length(ωs))
TOF_n2 = zeros(length(E_grid), length(ωs))
##

for (i, p) in enumerate(particles)
    TOF_n1[i, :] = WPI_TOF(ωs, p, plasma; n=1, wave_launch_time=wave_chirp)
    TOF_n2[i, :] = WPI_TOF(ωs, p, plasma; n=2, wave_launch_time=wave_chirp)
end


##
fig = Figure(size=(1200, 500))

Label(fig[0, 1:3],
    "Wave-particle interaction time-of-flight";
    fontsize=16, font=:bold
)

ax1 = Axis(fig[1, 1];
    xlabel = "Time [s]",
    ylabel = "E [keV]",
    title  = "n = 1"
)

ax2 = Axis(fig[1, 2];
    xlabel = "Time [s]",
    ylabel = "E [keV]",
    title  = "n = 2"
)

colors = [:blue, :red, :green, :orange, :purple]

for (i, ω) in enumerate(ωs)
    label = "ω = $(round(ω/plasma.Ω_e[1], digits=2))Ωe"
    lines!(ax1, TOF_n1[:, i], E_grid ./ 1e3; color=colors[i], label=label)
    lines!(ax2, TOF_n2[:, i], E_grid ./ 1e3; color=colors[i], label=label)
end

params = """
    Parameters:
    L = $L
    α = $(round(rad2deg(acos(abs(μ))), digits=1))°
    n_e = $(n_e0/1e6) cm⁻³
    """

Label(fig[1, 3],
    params;
    fontsize   = 10,
    tellwidth  = false,
    tellheight = false,
    halign     = :left,
    valign     = :bottom  # ← bottom of the same cell as legend
)

# Link axes so they share the same limits
linkyaxes!(ax1, ax2)
linkxaxes!(ax1, ax2)

Legend(fig[1, 3], ax1)
