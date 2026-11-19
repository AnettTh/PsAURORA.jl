using AURORA
using CairoMakie

## Define method
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    t;
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)

## Define magnetic fieldline
fieldline = StationFieldLine(magnetic_field, ref, t)

## Make figure
fig = Figure(size=(800, 400))
ax  = Axis(fig[1, 1];
    xlabel = "Rₚₑᵣₚ [RE]",
    ylabel = "Z (SM) [RE]",
    title  = "Electron density for $(station) at $t",
    aspect = DataAspect(),
)

## Colormap and global density range for consistent colorbar
log_n = log10.(n_e)
clims = extrema(log_n)

# Meridional projection, signed so the nightside is negative
xs = fieldline.xs ./ RE
ys = fieldline.ys ./ RE
zs = fieldline.zs ./ RE

ϕ_eq   = fieldline.ϕ_eq
r_sign = (π < ϕ_eq ≤ 2π) ? -1 : 1
r_perp = @. r_sign * sqrt(xs^2 + ys^2)

## Earth
θ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)

## Field line coloured by density (both hemispheres are already in the arrays)
lines!(ax, r_perp, zs; color=log_n, colormap=:turbo, colorrange=clims, linewidth=5)

Colorbar(fig[1, 2];
    colormap   = :turbo,
    limits     = clims,
    label      = "log₁₀(nₑ) [m⁻³]",
    tellheight = false,
)

xlims!(ax, 0, -10)
save("src/WPI/test_scripts/density/figures/11192026_TsyganenkoDentonDensity.png", fig)
