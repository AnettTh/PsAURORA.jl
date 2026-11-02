using AURORA
using JLD2
using CairoMakie

# [E, μ, ω, L, Ψ₀]
data = load("src/WPI/data/02112026_Psi0test_MLT4.0.jld2")

function nearest(grid, target)
    i = argmin(abs.(grid .- target))
    return grid[i], i
end

## Choose the fixed parameters by value
α_fix, iα = nearest(data["α_grid"], deg2rad(2))
L_fix, iL = nearest(data["L_grid"], 4.0)

xaxis   = data["ω_frac"]
yaxis   = data["E_grid"] ./ 1e3                 # [keV]
Ψ0_grid = data["Ψ0_grid"]
n_Ψ     = length(Ψ0_grid)

# Which (L, Ψ0) slices actually finished
done = data["done"][iL, :]
all(done) || @warn "Unfinished Ψ0 slices (they will be empty)" findall(!, done)

## Slice: [E, ω, Ψ0] at fixed α and L
TOF = data["TOF"][:, iα, :, iL, :]                # size (9, 10, 8)

panel_data(k)  = TOF[:, :, k]'                         # (n_ω × n_E), as heatmap! expects
panel_title(k) = "Ψ₀ = $(round(rad2deg(Ψ0_grid[k]), digits=0))°"

## Shared colour range (NaN and bounce zeros are excluded from the range only)
vals = [v for k in 1:n_Ψ for v in panel_data(k) if !isnan(v) && v != 0.0]
isempty(vals) && error("No resonant, non-bouncing particles in this slice")
colorrange = extrema(vals)

## Figure
n_rows, n_cols = 2, 4
fig = Figure(size=(1500, 700))
Label(fig[0, 1:n_cols],
    "Time-of-flight, α = $(round(rad2deg(α_fix), digits=1))°, L = $L_fix";
    fontsize=18, font=:bold, tellwidth=false)

panel_axes = Matrix{Axis}(undef, n_rows, n_cols)
hms        = Vector{Any}(undef, n_Ψ)

for k in 1:n_Ψ
    r, c = fldmod1(k, n_cols)
    ax = Axis(fig[r, c];
        title  = panel_title(k),
        xlabel = "ω/Ωₑ",
        ylabel = "E [keV]",
    )
    panel_axes[r, c] = ax

    hms[k] = heatmap!(ax, xaxis, yaxis, panel_data(k);
        colormap   = :turbo,
        colorrange = colorrange,
        nan_color  = :lightgray,    # no resonance
        lowclip    = :black,        # bounce zeros
    )
end

Colorbar(fig[1:n_rows, n_cols + 1], hms[1]; label="TOF [s]")

## Link axes and hide redundant labels
linkaxes!(vec(panel_axes)...)
for r in 1:n_rows, c in 1:n_cols
    r < n_rows && hidexdecorations!(panel_axes[r, c]; grid=false)
    c > 1      && hideydecorations!(panel_axes[r, c]; grid=false)
end

colgap!(fig.layout, 12)
rowgap!(fig.layout, 12)
fig
