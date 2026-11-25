using AURORA
using CairoMakie
using Dates


## Make the magnetic field and field line
station = :sodankyla
ref = LTPReference(station)
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)
fieldline = StationFieldLine(magnetic_field, ref, t)

density_model = DentonDensity(fieldline)
## Define variables for test-run
N = electron_density(density_model, fieldline.position_eq)

T_parallel = 20e3
A = 0.0         # Anisotropy
Δ = 0.0         # 0.0 means fully emptied losscone
r = 2.0         #
M = 20          # Number of subtractions

α_lc = losscone_angle(fieldline)

##
a_par = sqrt(2 * T_parallel * eV_in_J / mₑ)
v_par_grid  = range(-4 * a_par, 4 * a_par,  length=300)
v_perp_grid = range(-4 * a_par, 4 * a_par, length=300)


##
F = [subtracted_bimaxwellian(vpa, vpe, N, T_parallel, A, α_lc, Δ, r, M)
    for vpa in v_par_grid, vpe in v_perp_grid]


## Make plot
fig = Figure(size=(700, 600))
ax  = Axis(fig[1, 1];
    xlabel = "v∥ / c",
    ylabel = "v⊥ / c",
    title  = "Subtracted bi-Maxwellian \n
    T∥=$(round(T_parallel*1e-3, digits=1)) keV, A=$(round(A, digits=1)) keV, Δ=$Δ, r=$r",
    aspect=DataAspect()
)

## Meaningful range
F_log = log10.(max.(F, maximum(F) * 1e-20))
clims = (-30, maximum(F_log))
xlims!(ax, -1, 1)
ylims!(ax, -1, 1)
hm = heatmap!(ax,
    collect(v_par_grid) / c₀,
    collect(v_perp_grid) / c₀,
    F_log';
    colormap = :turbo,
    colorrange = clims
)

Colorbar(fig[1, 2], hm; label="log(Phase space density [m⁻⁶s³])", tellheight=false)
save("src/WPI/test_scripts/distributions/figures/SBM_Liu.png", fig)
