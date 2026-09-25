using AURORA
using CairoMakie

# Bin centers
L_centers = [3.5, 3.9, 4.4, 4.9, 5.5, 6.2, 7.0, 7.8]

# Bin edges (midpoints between centers, plus outer edges)
L_edges = [
    3.3,   # 3.5 - 0.2
    3.7,   # midpoint 3.5-3.9
    4.15,  # midpoint 3.9-4.4
    4.65,  # midpoint 4.4-4.9
    5.2,   # midpoint 4.9-5.5
    5.85,  # midpoint 5.5-6.2
    6.6,   # midpoint 6.2-7.0
    7.4,   # midpoint 7.0-7.8
    8.2,   # 7.8 + 0.4
]

# Table values
#              3.5,   3.9,   4.4,   4.9,  5.5,  6.2   7.0, 7.8
n_e0_vals = [530. , 380. , 230. , 140. , 83. , 39. , 15. , 7.7]  # in cm⁻³
α_vals    = [  0.2,   0.4,   0.8,   0.9,  0.8,  1.3,  2.1, 1.6]
L_α_vals  = [  8.1,   5.9,   4.8,   5.2,  6.4,  5.5,  4.8, 6.1]

##
function ne_denton(L, λ; SI::Bool=true)

    R = L .* RE .* cos.(λ).^2

    # Find the right bin
    i = searchsortedfirst(L_edges, L) - 1
    @show L, i, length(n_e0_vals)  # debug
    i = clamp(i, 1, length(n_e0_vals))

    # Convert to SI-units
    if SI
        n_e0 = n_e0_vals[i] * 1e6   # [m⁻³]
    else
        n_e0 = n_e0_vals[i]
    end

    α = α_vals[i]

    return n_e0 .* (L.*RE ./ R).^α
end


##
L_vals = [3, 4, 5, 6, 7, 8]
colors = [:blue, :red, :green, :orange, :purple, :brown]

fig = Figure()
ax  = Axis(fig[1, 1];
    xlabel = "R [RE]",
    ylabel = "nₑ [m⁻³]",
    title  = "Electron density vs radial distance",
    yscale = log10
)

for (L, color) in zip(L_vals, colors)
    # R ranges from Earth surface to beyond L*RE
    λ_max  = acos(sqrt((RE + z_ionosphere) / (L * RE)))
    λ_grid = range(0.0, λ_max, length=500)

    n_e = ne_denton.(L, λ_grid; SI=true)
    R_grid = @. L * RE * cos(λ_grid)^2

    lines!(ax, R_grid ./ RE, n_e; color=color, label="L = $L")
end

axislegend(ax, position=:rt)
save("src/WPI/diagnostic_figures/simple_denton_density.png", fig)


##
fig = Figure(size=(800, 700))
ax  = Axis(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Denton model for electron density in the magnetosphere",
    aspect = DataAspect()
)

# Earth
θ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)

# Colormap and global density range for consistent colorbar
n_e_all = Float64[]
for L in L_vals
    λ_max  = acos(sqrt((RE + z_ionosphere) / (L * RE)))
    λ_grid = range(0.0, λ_max, length=500)
    append!(n_e_all, log10.(ne_denton.(L, λ_grid; SI=true)))
end
clims = (minimum(n_e_all), maximum(n_e_all))

for L in L_vals
    λ_max  = acos(sqrt((RE + z_ionosphere) / (L * RE)))
    λ_grid = range(0.0, λ_max, length=500)

    n_e    = ne_denton.(L, λ_grid; SI=true)
    R_grid = @. L * RE * cos(λ_grid)^2
    x      = @. -R_grid * cos(λ_grid) / RE   # nightside → negative x
    z      = @.  R_grid * sin(λ_grid) / RE

    # Both hemispheres
    lines!(ax, x,  z; color=log10.(n_e), colormap=:turbo, colorrange=clims, linewidth=3)
    lines!(ax, x, -z; color=log10.(n_e), colormap=:turbo, colorrange=clims, linewidth=3)
end

Colorbar(fig[1, 2];
    colormap = :turbo,
    limits   = clims,
    label    = "log₁₀(nₑ) [m⁻³]",
    tellheight = false
)

xlims!(0, -8.5)
