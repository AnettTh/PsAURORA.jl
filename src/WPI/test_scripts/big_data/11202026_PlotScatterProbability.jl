using AURORA
using CairoMakie
using JLD2

@load "src/WPI/data/11-20-2026/CollectScatterProbability.jld2"

## Count how many of each energy and α₀ have precipitated
precipitated = @. isfinite(TOF) & (TOF > 0)
bounced = @. TOF == 0.0
resonant = @. !isnan(TOF)

collapse(A) = dropdims(sum(A; dims=(3, 4)); dims=(3, 4))

# Of size [E, α]
n_precip = collapse(precipitated)
n_bounced = collapse(bounced)
n_resonant = collapse(resonant)
n_total = size(TOF, 3) * size(TOF, 4)

probability = n_precip ./ n_total

function wilson(k, n; z=1.96)
    p      = k / n
    denom  = 1 + z^2 / n
    center = (p + z^2 / (2n)) / denom
    half   = z * sqrt(p * (1 - p) / n + z^2 / (4n^2)) / denom
    return center - half, center + half
end

lo = similar(probability); hi = similar(probability)
for I in CartesianIndices(probability)
    lo[I], hi[I] = wilson(n_precip[I], n_total)
end

## Make the figure
fig = Figure(size=(700, 110 * 6 + 100))
Label(
    fig[0, 1],
    text="Probability of precipitation for $(params.station) at $(params.time)",
    tellwidth=false,
    fontsize=18,
    font=:bold,
    )

E_keV = collect(E_grid) ./ 1e3
α_deg = collect(α0_grid)                     # must be in degrees
nE    = length(E_keV)

# Grid shape: as square as possible, one cell per energy
n_rows, n_cols = nE, 1
panel_axes = Matrix{Union{Axis, Nothing}}(nothing, n_rows, n_cols)

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

for iE in 1:nE
    ax = Axis(fig[iE, 1];
        xlabel = "α₀ [°]",
        #title  = "$(round(E_keV[iE], digits=1)) keV",
    )
    text!(ax, 0.98, 0.92;
        text     = "$(round(E_keV[iE], digits=1)) keV",
        space    = :relative,          # 0 to 1 across the axis, independent of the data limits
        align    = (:right, :top),
        fontsize = 14,
        font = :bold,
        color    = :black,
    )
    panel_axes[iE] = ax
    c = colors[mod1(iE, length(colors))]

    band!(ax, α_deg, lo[iE, :], hi[iE, :]; color=(c, 0.25))
    lines!(ax, α_deg, probability[iE, :]; color=colors[mod1(iE, length(colors))])
    ax.xticks = [0.1; 1:1:10]
end

# Same limits everywhere, and no repeated labels
linkaxes!(panel_axes...)
for iE in 1:nE-1
    hidexdecorations!(panel_axes[iE]; grid=false)
end
ylims!(panel_axes[1], 0, 0.25)
xlims!(panel_axes[1], 0.1, 10)

colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)
save("src/WPI/test_scripts/big_data/figures/11202026_ScatteringProbability.png", fig)
