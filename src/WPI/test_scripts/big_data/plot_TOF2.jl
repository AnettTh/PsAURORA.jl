using AURORA
using CairoMakie
using JLD2

#===

Make figures with one ax per E and α(ω), for given L

===#

## Ψ0, E_grid, L_grid, μ_grid, Bw, α_grid, iL_done, use_boris, ϕ_eq, ω_frac, TOF, Ω_e_eq
data = load("src/WPI/data/02112026_TOF_4D_MLT4.0.jld2")

size6x3 = (1300, 700)
n_rows, n_cols = 2, 3
n_panels = n_rows * n_cols

## Choosing my favorite values
function nearest(grid, target)
    i = argmin(abs.(grid .- target))
    return grid[i], i
end
favorite_L, i_favorite_L = nearest(data["L_grid"], 6.0)

## Gather the data
TOF = data["TOF"][:, :, :, i_favorite_L]         #[E, μ, ω, L]
x_axis = data["ω_frac"]
y_axis = data["μ_grid"]

## Find the panel indices, n_panels indices spread evenly over μ_grid
E_grid = data["E_grid"]
i_E_panels = round.(Int, range(1, length(E_grid), length=6))
E_panels = E_grid[i_E_panels]

## Functions to get data- and title to use on panel k
function panel_data(k)
    return view(TOF, i_E_panels[k], :, :)
end

function panel_title(k)
    return "E = $(round(E_grid[i_E_panels[k]] / 1e3, digits=1)) keV"
end

## Set common colorrange
all_vals = [v for k in 1:n_panels for v in panel_data(k)]
colorrange = extrema(all_vals)

## Initialize figure
fig2 = Figure(size=size6x3)
axes = Matrix{Axis}(undef, n_rows, n_cols)
tof = nothing

# Fill subfigures
for k in 1:n_panels
    r, c = fldmod1(k, n_cols)
    ax = Axis(fig2[r, c];
        title  = panel_title(k),
        xlabel = "ω/Ωₑ",
        ylabel = "α [°]",
    )
    axes[r, c] = ax

    # Transposed to get correct x- and y-axis
    tof = heatmap!(ax, x_axis, y_axis, panel_data(k)';
        colormap   = :turbo,
        colorrange = colorrange,
    )
end

# Make colorbar, link axes and hide redundant labels
Colorbar(fig2[1:n_rows, n_cols + 1], tof; label="TOF [s]")
used = filter(!isnothing, vec(axes))
linkaxes!(used...)
for r in 1:n_rows, c in 1:n_cols
    isnothing(axes[r, c]) && continue
    r < n_rows && !isnothing(axes[r + 1, c]) && hidexdecorations!(axes[r, c]; grid=false)
    c > 1 && hideydecorations!(axes[r, c]; grid=false)
end

colgap!(fig2.layout, 12)
rowgap!(fig2.layout, 12)
fig2
