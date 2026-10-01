using AURORA
using CairoMakie

##
L_vals = [3, 4, 5, 6, 7, 8]
ϕ = 0.0
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

    R_max = find_R_max(dipole_field, L, 0.0, ϕ)

    ne_func = (L, λ) -> denton_density_model(L, λ, ϕ, dipole_field, R_max)
    n_e     = ne_func.(L, λ_grid)

    #n_e = denton_density_model.(L, λ_grid, ϕ, dipole_field, R_max; SI=true)
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

    R_max = find_R_max(dipole_field, L, 0.0, ϕ)
    ne_func = (l, λ) -> denton_density_model(l, λ, ϕ, dipole_field, R_max)
    append!(n_e_all, log10.(ne_func.(L, λ_grid)))
    #append!(n_e_all, log10.(denton_density_model.(L, λ_grid, ϕ, dipole_field, R_max; SI=true)))
end
clims = (minimum(n_e_all), maximum(n_e_all))

for L in L_vals
    λ_max  = acos(sqrt((RE + z_ionosphere) / (L * RE)))
    λ_grid = range(0.0, λ_max, length=500)

    R_max = find_R_max(dipole_field, L, 0.0, ϕ)
    ne_func = (L, λ) -> denton_density_model(L, λ, ϕ, dipole_field, R_max)
    n_e     = ne_func.(L, λ_grid)
    #n_e    = denton_density_model.(L, λ_grid, ϕ, dipole_field, R_max; SI=true)
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
