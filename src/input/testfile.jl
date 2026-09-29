using AURORA
using CairoMakie

## Define variables for test-run
T_par = 1e2
T_perp = 1e4
N = 1e6
Δ = 0.5
β = 0.2
α_lc = deg2rad(3.0)

a_par  = sqrt(2 * T_par  * eV_in_J / mₑ)
a_perp = sqrt(2 * T_perp * eV_in_J / mₑ)

v_par_grid  = range(-5 * a_par,  5 * a_par,  length=300)
v_perp_grid = range(-5 * a_perp, 5 * a_perp, length=300)

F = [subtracted_bimaxwellian(vpar, vperp, N, a_par, a_perp, Δ, β, relativistic=false)
     for vperp in v_perp_grid, vpar in v_par_grid]

v_lc = collect(v_perp_grid) ./ tan(α_lc)


# TODO: Figure out why my perpendicular particles precipitate!!
## Make plot
fig = Figure(size=(650, 600))
ax  = Axis(fig[1, 1];
    xlabel = "v∥",
    ylabel = "v⊥",
    title  = "Loss cone bi-Maxwellian \n
    T∥=$(round(T_par*1e-3, digits=1)) keV, T⊥=$(round(T_perp*1e-3, digits=1)) keV, Δ=$Δ, β=$β",
    aspect=DataAspect()
)

## Meaningful range
clims = (6e-22, 5e-18)

hm = heatmap!(ax,
    collect(v_par_grid)  ./ a_par,
    collect(v_perp_grid) ./ a_perp,
    F;
    colormap = :turbo,
    colorrange = clims
)

#lines!(ax,  v_lc ./ a_par, collect(v_perp_grid) ./ a_perp;
#    color=:white, linestyle=:dash, linewidth=2)
#lines!(ax, -v_lc ./ a_par, collect(v_perp_grid) ./ a_perp;
#    color=:white, linestyle=:dash, linewidth=2)

Colorbar(fig[1, 2], hm; label="Phase space density [m⁻⁶s³]", tellheight=false)

##
xlims!(-5, 5)
ylims!(-5, 5)
