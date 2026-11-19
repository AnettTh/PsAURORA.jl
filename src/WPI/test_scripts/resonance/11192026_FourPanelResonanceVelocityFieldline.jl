using AURORA
using CairoMakie
using Dates

## Field, station line, particle and plasma
t       = DateTime("2020-12-21T02:00:00")
station = :sodankyla

magnetic_field = TsyganenkoMagneticField(t; pdyn=2.0, dst=-28.0, byimf=0.0, bzimf=2.0)
fieldline      = StationFieldLine(magnetic_field, LTPReference(station), t)

E_eV = 30e3
α0   = 0.5
μ0   = -cos(deg2rad(α0))
particle = ParticleState(E_eV, μ0, fieldline; relativistic=true)

density = DentonDensity(fieldline)
plasma  = PlasmaState(fieldline, density)

## Frequencies, one per panel
θ       = deg2rad(0.0)
n       = 1
ω_fracs = [0.1, 0.2, 0.3, 0.4]
Ω_eq    = minimum(plasma.Ω_e)

## Resonance velocity along the line for one frequency
function resonance_velocity(ω)
    return map(eachindex(plasma.Ω_e)) do i
        Ω_e  = plasma.Ω_e[i]
        ω_pe = plasma.ω_pe[i]

        k = try
            dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)
        catch
            # Outside valid range
            return NaN
        end
        k_par = k * cos(θ)
        iszero(k_par) && return NaN

        # Compute and remove faster-than-light results
        V = abs((ω - n * Ω_e / particle.γ) / k_par) / c₀
        return V > 1.0 ? NaN : V
    end
end

V_R = [resonance_velocity(f * Ω_eq) for f in ω_fracs]

## Geometry: meridional projection, signed so the nightside is negative
xs = fieldline.xs ./ RE
ys = fieldline.ys ./ RE
zs = fieldline.zs ./ RE

ϕ_eq   = fieldline.ϕ_eq
r_sign = (π < ϕ_eq ≤ 2π) ? -1 : 1
r_perp = @. r_sign * sqrt(xs^2 + ys^2)

## Shared colour range over all four panels
vals = [v for V in V_R for v in V if isfinite(v) && v > 0]
isempty(vals) && error("No physical resonance velocity for these settings")
clims = (minimum(vals), maximum(vals))

## Figure
n_rows, n_cols = 2, 2
fig = Figure(size=(900, 500))
panel_axes = Matrix{Axis}(undef, n_rows, n_cols)

Φ = range(0, 2π, length=100)

for (k, ω_frac) in enumerate(ω_fracs)
    r, c = fldmod1(k, n_cols)
    ax = Axis(fig[r, c];
        xlabel = "Rₚₑᵣₚ [RE]",
        ylabel = "Z (SM) [RE]",
        title  = "ω = $(ω_frac) Ωₑ",
        aspect = DataAspect(),
    )
    panel_axes[r, c] = ax

    poly!(ax, Point2f.(cos.(Φ), sin.(Φ)); color=:black, strokecolor=:black, strokewidth=1)
    lines!(ax, r_perp, zs;
        color=V_R[k], colormap=:turbo, colorrange=clims, nan_color=:lightgray, linewidth=6)

    xlims!(ax, 0, -10)
end

# Same limits in all panels, and no repeated labels
linkaxes!(vec(panel_axes)...)
for r in 1:n_rows, c in 1:n_cols
    r < n_rows && hidexdecorations!(panel_axes[r, c]; grid=false)
    c > 1      && hideydecorations!(panel_axes[r, c]; grid=false)
end
colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

Colorbar(fig[1:n_rows, n_cols + 1];
    colormap   = :turbo,
    limits     = clims,
    label      = "v∥ / c₀",
    tellheight = false,
)

## Parameters shared by all panels, in the title row
Label(fig[0, 1:n_cols + 1],
    "Resonance velocity for $(station) field line at $t   " *
    "(n = $n, θ = $(round(rad2deg(θ), digits=1))°)";
    fontsize=15, font=:bold, tellwidth=false)

fig
#save("src/WPI/test_scripts/resonance/figures/11192026_FourPanelResonanceVelocityFieldline.png", fig)
