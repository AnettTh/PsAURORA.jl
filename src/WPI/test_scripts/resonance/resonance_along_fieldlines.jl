using AURORA
using CairoMakie


## Shared setup
L_vals = [3, 4, 5, 6, 7]
n      = 1
θ      = range(0, 2π, length=100)

fig = Figure(size=(1400, 700))

## ── Panel 1: Fixed ω, varying v∥ (via α) ─────────────────────────────────────
ax1 = Axis(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Resonance latitude — fixed ω, varying α",
    aspect = DataAspect()
)

ω_fixed  = 0.35 * plasma.Ω_e[1]
α_vals   = [1.0, 2.0, 3.0, 5.0]   # pitch angles in degrees
colors_α = [:blue, :red, :green, :orange]

poly!(ax1, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)

for (α_deg, color) in zip(α_vals, colors_α)
    for L in L_vals
        plasma_L = PlasmaState(λ_grid, ne_denton, dipole_field, Float64(L))
        p        = ParticleState(30e3, -cos(deg2rad(α_deg)), [L*RE, 0.0, 0.0], dipole_field)
        λ_max    = acos(sqrt((RE + z_ionosphere) / (L * RE)))
        λ_grid_L = range(0.0, λ_max * 0.99, length=500)

        λ_res = resonance_latitude(ω_fixed, p, plasma_L; n=n)
        isnothing(λ_res) && continue
        # Field line
        R_fl = @. L * RE * cos.(λ_grid_L)^2
        x_fl = @. -R_fl *  cos.(λ_grid_L) / RE
        z_fl = @.  R_fl *  sin.(λ_grid_L) / RE
        lines!(ax1, x_fl, z_fl; color=:lightgray, linewidth=1)
        lines!(ax1, x_fl, -z_fl; color=:lightgray, linewidth=1)

        # Resonance point
        R_res = L * RE * cos(λ_res)^2
        x_res = -R_res * cos(λ_res) / RE
        z_res =  R_res * sin(λ_res) / RE
        scatter!(ax1, [x_res], [z_res]; color=color, markersize=8,
                 label= L==L_vals[1] ? "α=$(α_deg)°" : nothing)
        scatter!(ax1, [x_res], [-z_res]; color=color, markersize=8)
    end
end

axislegend(ax1, position=:rt)

## ── Panel 2: Fixed v∥ (via α), varying ω ─────────────────────────────────────
ax2 = Axis(fig[1, 2];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Resonance latitude — fixed α, varying ω",
    aspect = DataAspect()
)

α_fixed  = 3.0   # degrees
ω_fracs  = [0.2, 0.3, 0.4, 0.5]
colors_ω = [:blue, :red, :green, :orange]

poly!(ax2, Point2f.(cos.(θ), sin.(θ)); color=:lightblue, strokecolor=:black, strokewidth=1)

for (ω_frac, color) in zip(ω_fracs, colors_ω)
    for L in L_vals
        plasma_L = PlasmaState(λ_grid, ne_denton, dipole_field, Float64(L))
        p        = ParticleState(30e3, -cos(deg2rad(α_fixed)), [L*RE, 0.0, 0.0], dipole_field)
        ω        = ω_frac * plasma_L.Ω_e[1]
        λ_max    = acos(sqrt((RE + z_ionosphere) / (L * RE)))
        λ_grid_L = range(0.0, λ_max * 0.99, length=500)

        λ_res = resonance_latitude(ω, p, plasma_L; n=n)
        isnothing(λ_res) && continue

        R_fl = @. L * RE * cos(λ_grid_L)^2
        x_fl = @. -R_fl * cos(λ_grid_L) / RE
        z_fl = @.  R_fl * sin(λ_grid_L) / RE
        lines!(ax2, x_fl, z_fl; color=:lightgray, linewidth=1)
        lines!(ax2, x_fl, -z_fl; color=:lightgray, linewidth=1)

        R_res = L * RE * cos(λ_res)^2
        x_res = -R_res * cos(λ_res) / RE
        z_res =  R_res * sin(λ_res) / RE
        scatter!(ax2, [x_res], [z_res]; color=color, markersize=8,
                 label= L==L_vals[1] ? "ω=$(ω_frac)Ωe" : nothing)
        scatter!(ax2, [x_res], [-z_res]; color=color, markersize=8)
    end
end

axislegend(ax2, position=:rt)

#
save("src/PsA/plots/resonance_fieldlines.png", fig)
