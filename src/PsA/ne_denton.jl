using AURORA

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
n_e0_vals = [530., 380., 230., 140., 83., 39., 15., 7.7]  # in cm⁻³
α_vals    = [0.2,  0.4,  0.8,  0.9,  0.8, 1.3, 2.1, 1.6]
L_α_vals  = [8.1,  5.9,  4.8,  5.2,  6.4, 5.5, 4.8, 6.1]

##
function ne_denton(L, λ; SI::Bool=true)

    R = L .* RE .* cos.(λ).^2

    # Find the right bin
    i = searchsortedfirst(L_edges, L) - 1
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

    n_e = density_model.(L, λ_grid; SI=true)
    R_grid = @. L * RE * cos(λ_grid)^2

    lines!(ax, R_grid ./ RE, n_e; color=color, label="L = $L")
end

axislegend(ax, position=:rt)
##
save("src/PsA/density_vs_R.png", fig)
