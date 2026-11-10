using AURORA
using CairoMakie

## Parameters
L_vals = [5.0, 6.0, 7.0, 8.0]
MLTs   = [4.0, 6.0, 8.0]
ϕs     = deg2rad.(mod.((MLTs .- 12) .* 15, 360))

## Trace field lines and compute density
function trace_with_density(B_func, ne_func, x0, y0, z0; ds=RE*0.05, n_steps=1000)
    xs, ys, zs, n_es = [x0], [y0], [z0], Float64[]

    x, y, z = x0, y0, z0
    for _ in 1:n_steps
        B    = B_func(Cartesian(x, y, z))
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
function collect_fieldlines(adjusted_model, ϕ)
    magnetic_field = TsyganenkoMagneticField(ϕ)
    return map(L_vals) do L
        density = DentonDensity(magnetic_field, L, ϕ; adjusted_model=adjusted_model)

        x0 = L * RE * cos(ϕ)
        y0 = L * RE * sin(ϕ)
        z0 = 0.0

        xs_f, ys_f, zs_f, ne_f = trace_with_density(magnetic_field, density, x0, y0, z0; ds= RE*0.05)
        xs_b, ys_b, zs_b, ne_b = trace_with_density(magnetic_field, density, x0, y0, z0; ds=-RE*0.05)

        xs  = [reverse(xs_b);  xs_f[2:end]]
        ys  = [reverse(ys_b);  ys_f[2:end]]
        zs  = [reverse(zs_b);  zs_f[2:end]]
        n_es = [reverse(ne_b); ne_f[2:end]]
        return (; xs, ys, zs, n_es)
    end
end

# ADDED: rows = (adjusted, original), columns = MLT
models = (true, false)
fieldline_data = [collect_fieldlines(adj, ϕ) for adj in models, ϕ in ϕs]

## Colorrange
all_ne = vcat([log10.(d.n_es) for set in fieldline_data for d in set]...)
clims  = extrema(filter(isfinite, all_ne))

## Plot — x-z projection
fig = Figure(size=(900, 500))
θ = range(0, 2π, length=100)

# Store axes to link later
axs = Matrix{Axis}(undef, length(models), length(MLTs))

for row in eachindex(models), col in eachindex(MLTs)
    model_name = models[row] ? "Adjusted" : "Original"
    ax = Axis(fig[row, col];
        xlabel = "ρ [RE]",
        ylabel = "Z [RE]",
        title  = "$(model_name), $(MLTs[col]) MLT",
        aspect = DataAspect()
    )
    axs[row, col] = ax

    poly!(ax, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)

    for d in fieldline_data[row, col]
        # Project onto x-z plane using r_perp = sqrt(x²+y²) with sign from x
        r_perp = @. -sqrt(d.xs^2 + d.ys^2) / RE   # negative = nightside
        z_plot = d.zs ./ RE
        lines!(
            ax,
            r_perp,
            z_plot;
            color=log10.(d.n_es),
            colormap=:turbo,
            colorrange=clims,
            linewidth=3
        )
    end

    xlims!(ax, 0, -9)
    ylims!(ax, -4, 4)
end

# Link all axes, then hide redundant decorations
linkaxes!(vec(axs)...)
for row in eachindex(models), col in eachindex(MLTs)
    row < length(models) && hidexdecorations!(axs[row, col]; grid=false)
    col > 1              && hideydecorations!(axs[row, col]; grid=false)
end

colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

Colorbar(
    fig[1:2, length(MLTs) + 1];
    colormap   = :turbo,
    limits     = clims,
    label      = "log₁₀(nₑ) [m⁻³]",
    tellheight = false
)

Label(
    fig[0, 1:length(MLTs) + 1],
    "Electron density on Tsyganenko field lines";
    fontsize=18,
    font=:bold,
    tellwidth=false
)

##
save("src/WPI/test_scripts/density/figures/tsyg_density_fieldlines.png", fig)
