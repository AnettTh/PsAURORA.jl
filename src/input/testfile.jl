using AURORA
using CairoMakie

## Define variables for test-run
T_parallel = 10e3
N = 1e6
A = 1.0
Δ = 0.0
β = 0.8

a_par = sqrt(2 * T_parallel * eV_in_J / mₑ)
v_par_grid  = range(-4 * a_par, 4 * a_par,  length=300)
v_perp_grid = range(0, 4 * a_par, length=300)


##
F = [subtracted_bimaxwellian_Liu(vpar, vperp, N, T_parallel, A, Δ, β)
     for vperp in v_perp_grid, vpar in v_par_grid]


## Make plot
fig = Figure(size=(650, 600))
ax  = Axis(fig[1, 1];
    xlabel = "v∥",
    ylabel = "v⊥",
    title  = "Loss cone bi-Maxwellian \n
    T∥=$(round(T_parallel*1e-3, digits=1)) keV, A=$(round(A, digits=1)) keV, Δ=$Δ, β=$β",
    aspect=DataAspect()
)

## Meaningful range
F_log = log10.(max.(F, maximum(F) * 1e-6))
clims = (minimum(F_log), maximum(F_log))

hm = heatmap!(ax,
    collect(v_par_grid),
    collect(v_perp_grid),
    F_log;
    colormap = :turbo,
    colorrange = clims
)

Colorbar(fig[1, 2], hm; label="log(Phase space density [m⁻⁶s³])", tellheight=true)
