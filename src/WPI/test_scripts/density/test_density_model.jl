using AURORA
using CairoMakie

##
L_vals = [3]#, 4, 5, 6, 7, 8]
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
#magnetic_field = DipoleMagneticField()
colors = [:blue, :red, :green, :orange, :purple, :brown]


##TESTING
model = DentonDensity(
    magnetic_field,
    6.5,
    ϕ_eq;
    hemisphere=:north,
    adjusted_model=true,
)



##
fig = Figure()
ax  = Axis(
    fig[1, 1];
    xlabel = "R [RE]",
    ylabel = "nₑ [m⁻³]",
    title  = "Electron density vs radial distance",
    yscale = log10
)

##
for (L, color) in zip(L_vals, colors)
    # Using a dipole approximation, what is the maximum latitude reached?
    λ_max = acos(sqrt((RE + z_ionosphere) / (L * RE)))

    # Make grid to plot for
    λ_grid = range(0, λ_max, length=500)
    R_grid = @. L * RE * cos(λ_grid)^2

    # Define the density model
    model = DentonDensity(magnetic_field, Float64(L), ϕ_eq; adjusted_model=true)
    #n_e = model.(L, λ_grid)

    #lines!(ax, R_grid ./ RE, n_e; color=color, label="L = $L")
end

#axislegend(ax, position=:rt)
