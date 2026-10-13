using AURORA
using CairoMakie
using ProgressMeter

# TODO: figure out how to use this
## Define variables for test-run
T_parallel = 10e3
N = 1e6
A = 1.0
Δ = 0.2
r = 2.0
M = 20
α_lc = deg2rad(3)

## Make particle
L = 4
ϕ_eq = deg2rad(240.0)
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define grids
E_grid = range(1e3, 40e3, length=50)
α_grid = range(deg2rad(0.5), deg2rad(30), length=30)   # pitch angles to sample

## Make particles
p_relativistic = Matrix{ParticleState}(undef, length(E_grid), length(α_grid))

@showprogress for (j, α) in enumerate(α_grid)
    for (i, E) in enumerate(E_grid)
        p_relativistic[i, j] = ParticleState(
            E,
            -cos(α),
            r0,
            magnetic_field,
            ϕ_eq;
            relativistic=true,
            precomputed_magnetic_field=false
        )
    end
end

##
a_par = sqrt(2 * T_parallel * eV_in_J / mₑ)
v_par_grid  = range(-5 * a_par, 5 * a_par,  length=300)
v_perp_grid = range(-5 * a_par, 5 * a_par, length=300)

##
flux = zeros(length(E_grid))
@showprogress for (i, p) in enumerate(p_relativistic)
    w = distribution_weight(p, N, T_parallel, A, α_lc, Δ, r, M)
    w < threshold && continue

    tof = TOF_matrix_relativistic[i, :]   # TOF for each frequency

    # Weight contribution by phase space density
    flux[i] = w * mean(filter(!isnan, tof))
end

### Make plot
#fig = Figure(size=(700, 600))
#ax  = Axis(fig[1, 1];
#    xlabel = "v∥",
#    ylabel = "v⊥",
#    title  = "Subtracted bi-Maxwellian \n
#    T∥=$(round(T_parallel*1e-3, digits=1)) keV, A=$(round(A, digits=1)) keV, Δ=$Δ, r=$r",
#    aspect=DataAspect()
#)
#
### Meaningful range
#F_log = log10.(max.(F, maximum(F) * 1e-6))
#clims = (minimum(F_log), maximum(F_log))
#
#hm = heatmap!(ax,
#    collect(v_par_grid),
#    collect(v_perp_grid),
#    F_log';
#    colormap = :turbo,
#    colorrange = clims
#)
#hm_mirror = heatmap!(ax,
#    collect(v_par_grid),
#    -collect(v_perp_grid),
#    F_log';
#    colormap = :turbo,
#    colorrange = clims)

#Colorbar(fig[1, 2], hm; label="log(Phase space density [m⁻⁶s³])", tellheight=false)
