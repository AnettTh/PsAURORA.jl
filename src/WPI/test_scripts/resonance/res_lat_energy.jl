using AURORA
using CairoMakie

## Parameters
L_vals  = [4.0, 5.0, 6.0, 7.0, 8.0]
μ       = -cos(deg2rad(3.0))
ω_frac  = 0.35   # fraction of equatorial Ωe
n       = 1
E_grid  = range(1e3, 100e3, length=200)
colors  = [:blue, :red, :green, :orange, :purple]

## Compute
fig = Figure(size=(800, 600))
ax  = Axis(fig[1, 1];
    xlabel = "E [keV]",
    ylabel = "λ_res [°]",
    title  = "Resonance latitude vs energy",
    xscale = log10
)

for (L, color) in zip(L_vals, colors)
    r0       = [L*RE, 0.0, 0.0]
    λ_grid   = range(0.0, deg2rad(60), length=500)
    plasma_L = PlasmaState(λ_grid, ne_denton, tsyganenko_field, L)
    ω        = ω_frac * plasma_L.Ω_e[1]

    λ_res_grid = map(E_grid) do E
        p     = ParticleState(E, μ, r0, tsyganenko_field; relativistic=true)
        λ_res = resonance_latitude(ω, p, plasma_L; n=n)
        isnan(λ_res) ? NaN : rad2deg(λ_res)
    end

    lines!(ax, collect(E_grid)./1e3, λ_res_grid;
        color=color, label="L = $L")
end

# Add parameter text
text!(ax, 0.02, 0.98;
    text      = "ω = $(ω_frac)Ωe, α = $(round(rad2deg(acos(abs(μ))), digits=1))°, n=$n",
    space     = :relative,
    align     = (:left, :top),
    fontsize  = 11
)

axislegend(ax, position=:rb)

save("src/WPI/test_scripts/resonance/resonance_lat_vs_energy_log.png", fig)
