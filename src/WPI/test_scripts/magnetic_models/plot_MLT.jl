using AURORA
using CairoMakie
using LinearAlgebra
## Grid parameters
n_points = 50
L_vals   = [4.0, 5.0, 6.0, 7.0, 8.0]
φ_MLT    = deg2rad(120.0)   # 4 MLT

## Trace field lines
function trace_fieldline(B_func, x0, y0, z0; ds=RE*0.05, n_steps=1000)
    xs, ys, zs = [x0], [y0], [z0]
    Bs = Float64[]

    x, y, z = x0, y0, z0
    for _ in 1:n_steps
        B    = B_func(x, y, z)
        Bmag = norm(B)
        push!(Bs, Bmag)
        Bx, By, Bz = B ./ Bmag
        x += ds * Bx
        y += ds * By
        z += ds * Bz
        push!(xs, x); push!(ys, y); push!(zs, z)
        norm([x, y, z]) < RE        && break
        norm([x, y, z]) > 10*RE     && break
    end
    push!(Bs, Bs[end])
    return xs, ys, zs, Bs
end

## Collect field line data — seed at equator for each L
fieldline_data = map(L_vals) do L
    x0 = L * RE * cos(φ_MLT)
    y0 = L * RE * sin(φ_MLT)
    z0 = 0.0

    xs_f, ys_f, zs_f, Bs_f = trace_fieldline(tsyganenko_field, x0, y0, z0; ds= RE*0.05)
    xs_b, ys_b, zs_b, Bs_b = trace_fieldline(tsyganenko_field, x0, y0, z0; ds=-RE*0.05)

    xs = [reverse(xs_b); xs_f[2:end]]
    ys = [reverse(ys_b); ys_f[2:end]]
    zs = [reverse(zs_b); zs_f[2:end]]
    Bs = [reverse(Bs_b); Bs_f[2:end]]
    return (; xs, ys, zs, Bs)
end

## Colorrange
all_Bs = vcat([log10.(d.Bs) for d in fieldline_data]...)
clims  = extrema(filter(isfinite, all_Bs))

## Figure — two panels: x-z plane and x-y plane
fig = Figure(size=(1200, 600))

# Panel 1: x-z plane (meridional)
ax1 = Axis(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Tsyganenko field — meridional plane (4 MLT)",
    aspect = DataAspect()
)

ϕ = range(0, 2π, length=100)
poly!(ax1, Point2f.(cos.(ϕ), sin.(ϕ)); color=:lightblue, strokecolor=:black, strokewidth=1)

for d in fieldline_data
    lines!(ax1, d.xs./RE, d.zs./RE;
        color=log10.(d.Bs), colormap=:plasma, colorrange=clims, linewidth=2)
end

# Panel 2: x-y plane (equatorial)
ax2 = Axis(fig[1, 2];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    title  = "Tsyganenko field — equatorial plane (4 MLT)",
    aspect = DataAspect()
)

poly!(ax2, Point2f.(cos.(ϕ), sin.(ϕ)); color=:lightblue, strokecolor=:black, strokewidth=1)

for (L, d) in zip(L_vals, fieldline_data)
    lines!(ax2, d.xs./RE, d.ys./RE;
        color=log10.(d.Bs), colormap=:plasma, colorrange=clims, linewidth=2)
    # Mark seed point (equator)
    scatter!(ax2, [L*cos(φ_MLT)], [L*sin(φ_MLT)]; color=:white, markersize=8)
end

# Add MLT direction indicator
lines!(ax2, [0, 9*cos(φ_MLT)], [0, 9*sin(φ_MLT)];
    color=:white, linestyle=:dash, linewidth=1)
text!(ax2, 9.2*cos(φ_MLT), 9.2*sin(φ_MLT); text="4 MLT", color=:white)

Colorbar(fig[1, 3];
    colormap   = :plasma,
    limits     = clims,
    label      = "log₁₀(|B|) [T]",
    tellheight = false
)
##
save("src/WPI/diagnostic_figures/tsyganenko_fieldlines.png", fig)
