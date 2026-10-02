using AURORA
using CairoMakie

## Define variables for test-run
T_parallel = 10e3
N = 1e6
A = 1.0
Δ = 0.5
β = 0.8
r = 2.0
M = 20
α_lc = deg2rad(3)

a_par = sqrt(2 * T_parallel * eV_in_J / mₑ)
v_par_grid  = range(-5 * a_par, 5 * a_par,  length=300)
v_perp_grid = range(-5 * a_par, 5 * a_par, length=300)


##
#F = [subtracted_bimaxwellian_Liu(vpar, vperp, N, T_parallel, A, Δ, β)
#     for vperp in v_perp_grid, vpar in v_par_grid]
F = [subtracted_bimaxwellian(vpa, vpe, N, T_parallel, A, α_lc, Δ, r, M)
    for vpa in v_par_grid, vpe in v_perp_grid]


## Make plot
fig = Figure(size=(700, 600))
ax  = Axis(fig[1, 1];
    xlabel = "v∥",
    ylabel = "v⊥",
    title  = "Subtracted bi-Maxwellian \n
    T∥=$(round(T_parallel*1e-3, digits=1)) keV, A=$(round(A, digits=1)) keV, Δ=$Δ, β=$β",
    aspect=DataAspect()
)

## Meaningful range
F_log = log10.(max.(F, maximum(F) * 1e-6))
clims = (minimum(F_log), maximum(F_log))

hm = heatmap!(ax,
    collect(v_par_grid),
    collect(v_perp_grid),
    F_log';
    colormap = :turbo,
    colorrange = clims
)
#hm_mirror = heatmap!(ax,
#    collect(v_par_grid),
#    -collect(v_perp_grid),
#    F_log';
#    colormap = :turbo,
#    colorrange = clims)

Colorbar(fig[1, 2], hm; label="log(Phase space density [m⁻⁶s³])", tellheight=false)
