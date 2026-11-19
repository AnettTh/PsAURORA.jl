using AURORA
using CairoMakie
using Dates

## Field and station line
t       = DateTime("2020-12-21T02:00:00")
station = :sodankyla

magnetic_field = TsyganenkoMagneticField(t; pdyn=2.0, dst=-28.0, byimf=0.0, bzimf=2.0)
fieldline = StationFieldLine(magnetic_field, LTPReference(station), t)

xs, ys, zs = fieldline.xs ./ RE, fieldline.ys ./ RE, fieldline.zs ./ RE
logB = log10.(fieldline.Bs)
clims = extrema(filter(isfinite, logB))

# Distance along the meridian through the magnetic equator, and the MLT of that meridian
ϕ_eq = fieldline.ϕ_eq
ρ = @. xs * cos(ϕ_eq) + ys * sin(ϕ_eq)
mlt = longitude_to_MLT(ϕ_eq; degrees=true)
i_eq = fieldline.i_eq

## Make figure in both meridional plane and equatorial plane
fig = Figure(size=(800, 500))
θ   = range(0, 2π, length=100)

# Panel 1: meridional plane through the magnetic equator
ax1 = Axis(fig[1, 1];
    xlabel = "Distance along the equator meridian [RE]",
    ylabel = "Z (SM) [RE]",
    title  = "Meridional plane ($(round(mlt, digits=1)) MLT at the equator)",
    aspect = DataAspect(),
)
poly!(ax1, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)
lines!(ax1, ρ, zs; color=logB, colormap=:turbo, colorrange=clims, linewidth=3)

# Panel 2: equatorial plane (SM x-y)
ax2 = Axis(fig[1, 2];
    xlabel = "X (SM) [RE]",
    ylabel = "Y (SM) [RE]",
    title  = "Equatorial plane",
    aspect = DataAspect(),
)
poly!(ax2, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)
lines!(ax2, xs, ys; color=logB, colormap=:turbo, colorrange=clims, linewidth=3)

Colorbar(fig[1, 3];
    colormap   = :turbo,
    limits     = clims,
    label      = "log₁₀(|B|) [T]",
    tellheight = false,
)

Label(fig[0, 1:3], "Field line through $(station) at $(t)";
    fontsize=18, font=:bold, tellwidth=false)

xlims!(ax2, 0, -5)
xlims!(ax1, 0, 10)
save("src/WPI/test_scripts/magnetic_models/figures/11192026_TsyganenkoFieldLines.png", fig)
