using AURORA
using AURORA; c₀
using CairoMakie

## Define configurations to compare
configs = [
    (
        label  = "L=6.5, E=30keV, α=3°",
        particle = ParticleState(3e4, -cos(deg2rad(3)), [6.5*RE, 0.0, 0.0], dipole_field; relativistic=true),
        plasma   = PlasmaState(range(0.0, deg2rad(50), length=500), 1.8e7, dipole_field, 6.5)
    ),
    (
        label  = "L=6.5, E=30keV, α=3°",
        particle = ParticleState(3e4, -cos(deg2rad(3)), [6.5*RE, 0.0, 0.0], dipole_field; relativistic=true),
        plasma   = PlasmaState(range(0.0, deg2rad(50), length=500), 1.8e7, dipole_field, 6.5)
    ),
    (
        label  = "L=6.5, E=30keV, α=3°",
        particle = ParticleState(3e4, -cos(deg2rad(3)), [6.5*RE, 0.0, 0.0], dipole_field; relativistic=true),
        plasma   = PlasmaState(range(0.0, deg2rad(50), length=500), 1.8e7, dipole_field, 6.5)
    ),
]

## Plot
fig = Figure()
ax  = Axis(fig[1, 1], xlabel="ω [rad/s]", ylabel="λᵣₑₛ [deg]")

colors = [:blue, :red, :green, :orange, :purple]

for (i, cfg) in enumerate(configs)
    ω_grid  = range(0.25 * cfg.plasma.Ω_e[1], 0.5 * cfg.plasma.Ω_e[1], length=500)
    λ_res   = resonance_latitude(ω_grid, cfg.particle, cfg.plasma)
    λ_plot  = [isnothing(x) ? NaN : rad2deg(x) for x in λ_res]

    lines!(ax, ω_grid, λ_plot; color=colors[i], label=cfg.label)
end

axislegend(ax, position=:rt)
fig



## Define configurations to compare
configs = [
    (
        label  = "ω=0.4Ωe, α=3°, L=6.5" * "nₑ=1.5e7",
        ω_frac = 0.4,
        α_deg  = 1.0,
        L      = 6.5,
        n_e0   = 1.5e7,
    ),
    (
        label  = "ω=0.4Ωe, α=3°, L=6.5" * "nₑ=1.8e7",
        ω_frac = 0.4,
        α_deg  = 3.0,
        L      = 6.5,
        n_e0   = 1.8e7,
    ),
    (
        label  = "ω=0.4Ωe, α=3°, L=6.5" * "nₑ=2.1e7",
        ω_frac = 0.4,
        α_deg  = 5.0,
        L      = 6.5,
        n_e0   = 2.1e7,
    ),
]

E_grid = range(1e3, 100e3, length=200)   # [eV]

## Plot
fig = Figure()
ax  = Axis(fig[1, 1],
    xlabel = "E [keV]",
    ylabel = "λᵣₑₛ [deg]",
    xscale = log10
)

colors = [:blue, :red, :green, :orange, :purple]

for (i, cfg) in enumerate(configs)
    plasma   = PlasmaState(range(0.0, deg2rad(50), length=500), cfg.n_e0, dipole_field, cfg.L)
    ω        = cfg.ω_frac * plasma.Ω_e[1]
    μ        = -cos(deg2rad(cfg.α_deg))
    r0       = [cfg.L * RE, 0.0, 0.0]

    λ_res_grid = map(E_grid) do E
        p     = ParticleState(E, μ, r0, dipole_field; relativistic=true)
        λ_res = resonance_latitude(ω, p, plasma)
        isnothing(λ_res) ? NaN : rad2deg(λ_res)
    end

    lines!(ax, E_grid ./ 1e3, λ_res_grid; color=colors[i], label=cfg.label)
end

axislegend(ax, position=:rb)
fig
