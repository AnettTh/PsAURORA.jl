using AURORA
using CairoMakie
using Dates

## Field, station line, particle and plasma
t       = DateTime("2020-12-21T02:00:00")
station = :sodankyla

magnetic_field = TsyganenkoMagneticField(t; pdyn=2.0, dst=-28.0, byimf=0.0, bzimf=2.0)
fieldline      = StationFieldLine(magnetic_field, LTPReference(station), t)

E_eV = 30e3
α0 = 3.0
μ0 = -cos(deg2rad(α0))
particle = ParticleState(E_eV, μ0, fieldline; relativistic=true)

density = DentonDensity(fieldline)
plasma  = PlasmaState(fieldline, density)

## Resonance velocity along the line
θ      = 0.0
n      = 1
ω_frac = 0.2
Ω_eq   = minimum(plasma.Ω_e)
ω      = ω_frac * Ω_eq

V_R = map(eachindex(plasma.Ω_e)) do i
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
## Geometry: meridional projection, signed so the nightside is negative
xs = fieldline.xs ./ RE
ys = fieldline.ys ./ RE
zs = fieldline.zs ./ RE

ϕ_eq   = fieldline.ϕ_eq
r_sign = (π < ϕ_eq ≤ 2π) ? -1 : 1
r_perp = @. r_sign * sqrt(xs^2 + ys^2)

## Colour range
vals = filter(v -> isfinite(v) && v > 0, V_R)
clims = (minimum(vals), maximum(vals))

## Figure
fig = Figure(size=(700, 400))
ax  = Axis(fig[1, 1];
    xlabel = "Rₚₑᵣₚ [RE]",
    ylabel = "Z (SM) [RE]",
    title  = "Resonance velocity for $(station) field line at $t",
    aspect = DataAspect(),
)

Φ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(Φ), sin.(Φ)); color=:black, strokecolor=:black, strokewidth=1)

lines!(ax, r_perp, zs;
    color=V_R, colormap=:turbo, colorrange=clims, nan_color=:lightgray, linewidth=6)

## Parameter text
parameters_text = """
    n = $n
    ω = $(ω_frac) Ωₑ
    θ = $(round(rad2deg(θ), digits=1))°
    """

text!(ax, 0.02, 0.98;
    text     = parameters_text,
    space    = :relative,
    align    = (:left, :top),
    fontsize = 11,
)


Colorbar(fig[1, 2];
    colormap   = :turbo,
    limits     = clims,
    label      = "v∥ / c₀",
    tellheight = false,
)

xlims!(ax, 0, -10)
save("src/WPI/test_scripts/resonance/figures/11192026_ResonanceVelocityFieldline.png", fig)
