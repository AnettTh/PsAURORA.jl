using AURORA
using CairoMakie
using JLD2

#===

Make figures with one ax per α and E(ω), for given L

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
#favorite_μ, i_favorite_μ = nearest(data["μ_grid"], -cos(deg2rad(2)))
#favorite_E, i_favorite_E = nearest(data["E_grid"], 25*1e3)
#favorite_ω, i_favorite_ω = nearest(data["ω_frac"], 0.3)


## Gather the data
TOF = data["TOF"][:, :, :, i_favorite_L]         #[E, μ, ω, L]
x_axis = data["ω_frac"]
y_axis = data["E_grid"] ./ 1e3

## Find the panel indices, n_panels indices spread evenly over μ_grid
μ_grid = data["μ_grid"]
i_μ_panels = 1:1:n_panels # the six first
μ_panels = μ_grid[i_μ_panels]

## Functions to get data- and title to use on panel k
function panel_data(k)
    return view(TOF, :, i_μ_panels[k], :)
end

function panel_title(k)
    α = rad2deg(acos(abs(μ_grid[i_μ_panels[k]])))
    return "α = $(round(α, digits=1))°"
end

## Set common colorrange
all_vals = [v for k in 1:n_panels for v in panel_data(k)]
colorrange = extrema(all_vals)

## Initialize figure
fig1 = Figure(size=size6x3)
axes = Matrix{Axis}(undef, n_rows, n_cols)
tof = nothing

# Fill subfigures
for k in 1:n_panels
    r, c = fldmod1(k, n_cols)
    ax = Axis(fig1[r, c];
        title  = panel_title(k),
        xlabel = "ω/Ωₑ",
        ylabel = "E [keV]",
    )
    axes[r, c] = ax

    # Transposed to get correct x- and y-axis
    tof = heatmap!(ax, x_axis, y_axis, panel_data(k)';
        colormap   = :turbo,
        colorrange = colorrange,
    )
end

# Make colorbar, link axes and hide redundant labels
Colorbar(fig1[1:n_rows, n_cols + 1], tof; label="TOF [s]")
used = filter(!isnothing, vec(axes))
linkaxes!(used...)
for r in 1:n_rows, c in 1:n_cols
    isnothing(axes[r, c]) && continue
    r < n_rows && !isnothing(axes[r + 1, c]) && hidexdecorations!(axes[r, c]; grid=false)
    c > 1 && hideydecorations!(axes[r, c]; grid=false)
end

colgap!(fig1.layout, 12)
rowgap!(fig1.layout, 12)
fig1
