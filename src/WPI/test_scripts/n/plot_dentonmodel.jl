using AURORA
using CairoMakie

##
L_vals = [3, 4, 5, 6, 7, 8]
ϕ = deg2rad(120.0)
magnetic_field = TsyganenkoMagneticField(ϕ)
colors = [:blue, :red, :green, :orange, :purple, :brown]

fig = Figure()
ax  = Axis(fig[1, 1];
    xlabel = "R [RE]",
    ylabel = "nₑ [m⁻³]",
    title  = "Electron density vs radial distance for Tsyganenko field",
    yscale = log10
)

for (L, color) in zip(L_vals, colors)
    λ_max  = acos(sqrt((RE + z_ionosphere) / (L * RE)))
    λ_grid = range(0.0, λ_max, length=500)

    model  = DentonDensity(magnetic_field, Float64(L), ϕ)
    n_e    = model.(L, λ_grid)

    R_grid = @. L * RE * cos(λ_grid)^2

    lines!(ax, R_grid ./ RE, n_e; color=color, label="L = $L")
end

axislegend(ax, position=:rt)
save("src/WPI/diagnostic_figures/simple_denton_density_tsyg.png", fig)


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
    model  = DentonDensity(magnetic_field, Float64(L), ϕ)
    append!(n_e_all, log10.(model.(L, λ_grid)))
end
clims = (minimum(n_e_all), maximum(n_e_all))

for L in L_vals
    λ_max  = acos(sqrt((RE + z_ionosphere) / (L * RE)))
    λ_grid = range(0.0, λ_max, length=500)
    model  = DentonDensity(magnetic_field, L, ϕ)
    n_e    = model.(L, λ_grid)
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
