using AURORA
using CairoMakie

## Define configurations to compare
function make_configs(;
    E_eV   = 30e3,
    α_deg  = 3.0,
    L      = 6.5,
    n_e0   = 1.8e7,
)
    # Detect which parameter is a collection
    params = (; E_eV, α_deg, L, n_e0)
    varying = [(k, v) for (k, v) in pairs(params) if v isa AbstractArray]

    # If none vary, return a single config
    if isempty(varying)
        return [_make_config(E_eV, α_deg, L, n_e0)]
    end

    key, values = varying[1]
    return map(values) do val
        E    = key == :E_eV  ? val : E_eV
        α    = key == :α_deg ? val : α_deg
        l    = key == :L     ? val : L
        ne   = key == :n_e0  ? val : n_e0
        _make_config(E, α, l, ne)
    end
end

function _make_config(E_eV, α_deg, L, n_e0)
    r0 = [L*RE, 0.0, 0.0]
    return (
        label    = "L=$L, E=$(E_eV/1e3)keV, α=$(α_deg)°, n_e=$(n_e0/1e6)cm⁻³",
        particle = ParticleState(E_eV, -cos(deg2rad(α_deg)), r0, dipole_field),
        plasma   = PlasmaState(range(0.0, deg2rad(50), length=500), n_e0, dipole_field, L)
    )
end

##
configs = make_configs(L=range(3.5, 6.5, length=10))

## resonance latitude vs frequency
ω_grid = range(0.25 * configs[1].plasma.Ω_e[1], 0.5 * configs[1].plasma.Ω_e[1], length=500)

fig = Figure()
ax  = Axis(fig[1, 1];
    xlabel = "ω/Ωₑ",
    ylabel = "λᵣₑₛ [°]",
    title  = "Resonance latitude vs frequency"
)

colors = cgrad(:turbo, length(configs), categorical=true)

for (i, cfg) in enumerate(configs)
    ω_grid_i = range(
        0.25 * cfg.plasma.Ω_e[1],
        0.5  * cfg.plasma.Ω_e[1],
        length = 500
    )
    λ_res = resonance_latitude(ω_grid_i, cfg.particle, cfg.plasma)
    ω_norm = ω_grid_i ./ cfg.plasma.Ω_e[1]  # normalize by local Ωe

    lines!(ax, ω_norm, rad2deg.(λ_res);
        color = colors[i],
        label = "L = $(round(cfg.particle.L, digits=1))"
    )
end

Legend(fig[1, 2], ax)
save("src/WPI/test_scripts/L/res_vs_freq.png",fig)
