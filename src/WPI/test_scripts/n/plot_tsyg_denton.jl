using AURORA
using CairoMakie

## Parameters
L_vals = [4, 5, 6, 7, 8]
ϕ      = deg2rad(120.0)   # 4 MLT

## Trace field lines and compute density
function trace_with_density(B_func, ne_func, x0, y0, z0; ds=RE*0.05, n_steps=1000)
    xs, ys, zs, n_es = [x0], [y0], [z0], Float64[]

    x, y, z = x0, y0, z0
    for _ in 1:n_steps
        B    = B_func(x, y, z)
        Bmag = norm(B)
        iszero(Bmag) && break
        bx, by, bz = B ./ Bmag

        # Compute density at current position
        r   = sqrt(x^2 + y^2 + z^2)
        λ   = asin(z / r)
        L_i = r / (RE * cos(λ)^2)
        push!(n_es, ne_func(L_i, λ))

        x += ds * bx
        y += ds * by
        z += ds * bz

        push!(xs, x); push!(ys, y); push!(zs, z)
        norm([x, y, z]) < RE       && break
        norm([x, y, z]) > 10 * RE  && break
    end
    push!(n_es, n_es[end])  # match length
    return xs, ys, zs, n_es
end

## Collect data
fieldline_data = map(L_vals) do L
    R_max   = find_R_max(tsyganenko_field, Float64(L), 0.0, ϕ)
    ne_func = (L_i, λ) -> denton_density_model(L_i, λ, ϕ, tsyganenko_field, R_max)

    x0 = L * RE * cos(ϕ)
    y0 = L * RE * sin(ϕ)
    z0 = 0.0

    xs_f, ys_f, zs_f, ne_f = trace_with_density(tsyganenko_field, ne_func, x0, y0, z0; ds= RE*0.05)
    xs_b, ys_b, zs_b, ne_b = trace_with_density(tsyganenko_field, ne_func, x0, y0, z0; ds=-RE*0.05)

    xs  = [reverse(xs_b);  xs_f[2:end]]
    ys  = [reverse(ys_b);  ys_f[2:end]]
    zs  = [reverse(zs_b);  zs_f[2:end]]
    n_es = [reverse(ne_b); ne_f[2:end]]
    return (; xs, ys, zs, n_es)
end

## Colorrange
all_ne = vcat([log10.(d.n_es) for d in fieldline_data]...)
clims  = extrema(filter(isfinite, all_ne))

## Plot — x-z projection
fig = Figure(size=(800, 600))
ax  = Axis(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Electron density — Tsyganenko field lines (4 MLT)",
    aspect = DataAspect()
)

θ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)

for d in fieldline_data
    # Project onto x-z plane using r_perp = sqrt(x²+y²) with sign from x
    r_perp = @. -sqrt(d.xs^2 + d.ys^2) / RE   # negative = nightside
    z_plot = d.zs ./ RE
    lines!(ax, r_perp,  z_plot; color=log10.(d.n_es), colormap=:turbo, colorrange=clims, linewidth=3)
    #lines!(ax, r_perp, -z_plot; color=log10.(d.n_es), colormap=:turbo, colorrange=clims, linewidth=3)
end

Colorbar(fig[1, 2];
    colormap   = :turbo,
    limits     = clims,
    label      = "log₁₀(nₑ) [m⁻³]",
    tellheight = false
)

xlims!(ax, 0, -9)

##
save("src/WPI/diagnostic_figures/tsyganenko_density_fieldlines.png", fig)
