using AURORA
using CairoMakie
using JLD2
using Dates
using ProgressMeter

## Define configuration
name = "src/WPI/data/11-09-2026/footpoint_and_phase.jld2"
use_boris = true
relativistic = true
WPI = true
field_dependent = true
n = 1
adjusted_model = false
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
ref = LTPReference(:sodankyla)
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    ref;
    time=t,
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)

## Define the particle
α_grid = deg2rad.([0.5, 1.0, 1.5, 2.0])
μ_grid = -cos.(α_grid)

# Find the L-shell and equatorial position of the fieldline belonging to the station at the
# requested time
trace_L_shell = find_L_shell(magnetic_field, ref, t)
L = trace_L_shell.L
position_eq = trace_L_shell.position_eq

ϕ_eq = magnetic_field.ϕ
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))
E_grid = range(7e3, 40e3, length=10)

particles = Matrix{ParticleState}(undef, length(E_grid), length(μ_grid))
@showprogress for iα in eachindex(α_grid), (iE, E) in enumerate(E_grid)
    particles[iE, iα] = ParticleState(
        E,
        μ_grid[iα],
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=relativistic,
    )
end

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
lb_low, lb_high = 0.1, 0.5
density = DentonDensity(
    magnetic_field,
    8.2,
    ϕ_eq;
    duct_enhancement=duct_enhancement,
    adjusted_model=adjusted_model,
)

plasma = PlasmaState(λ_grid, density; lb_low=lb_low, lb_high=lb_high)

## Define the wave
ω_frac_used = 0.3
ωs = [ω_frac_used * minimum(plasma.Ω_e)]
amplitude = ConstantAmplitude(100e-12)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

Ψ0_grid = collect(range(0.0, 2π, length=5))[1:end-1]
waves = [
    WPI ? WaveState(modulation, amplitude, plasma; Ψ0=Ψ0, wave_direction=-1) : nothing
    for Ψ0 in Ψ0_grid]


## Calculate TOF and footpoints
TOF = fill(NaN, length(E_grid), length(μ_grid), length(Ψ0_grid))
FP  = fill(NaN, 3, length(E_grid), length(μ_grid), length(Ψ0_grid))

pairs_Eα = CartesianIndices((length(E_grid), length(μ_grid)))
progress = Progress(length(pairs_Eα))
Threads.@threads for k in eachindex(pairs_Eα)
    iE, iα = Tuple(pairs_Eα[k])

    for iΨ in eachindex(Ψ0_grid)
        result = WPI_TOF(
            ωs,
            particles[iE, iα],
            plasma;
            wave=waves[iΨ],
            use_boris=use_boris,
            field_dependent=field_dependent,
            n=n,
            r_start=position_eq,
            get_footpoint=true,
        )
        TOF[iE, iα, iΨ] = result.tof[1]
        FP[:, iE, iα, iΨ] = result.footpoint[:, 1]
    end

    next!(progress)
end


## Convert to East-North-Up coordinates with chosen station as reference
reference = LTPReference(station)
fp_enu = fill(NaN, size(FP))
for iΨ in axes(FP, 4), iα in axes(FP, 3), iE in axes(FP, 2)
    fp_enu[:, iE, iα, iΨ] .= FP_to_ENU(FP[:, iE, iα, iΨ], reference, t)
end


##
E_keV   = collect(E_grid) ./ 1e3

east  = fp_enu[1, :, :, :] ./ 1e3                      # CHANGED: n_E x n_α x n_Ψ, in km
north = fp_enu[2, :, :, :] ./ 1e3

## Figure: one axis per pitch angle, stacked
n_rows, n_cols = 2, 2                                  # CHANGED: grid shape
fig = Figure(size = (800, 400))                       # CHANGED: wider, shorter
axs = Matrix{Axis}(undef, n_rows, n_cols)              # CHANGED: matrix instead of vector

sc = nothing
for iα in eachindex(α_grid)
    r, c = fldmod1(iα, n_cols)                         # CHANGED: loop over pitch angles
    ax = Axis(fig[r, c];
        xlabel = "East [km]",
        ylabel = "North [km]",
        title  = "α = $(round(rad2deg(α_grid[iα]), digits=1)) deg",
        aspect = DataAspect(),
    )
    axs[r, c] = ax

    # One scatter per phase, coloured by energy
    for iΨ in eachindex(Ψ0_grid)
        sc = scatter!(ax, east[:, iα, iΨ], north[:, iα, iΨ];
            color      = E_keV,
            colormap   = :turbo,
            colorrange = extrema(E_keV),
            markersize = 12,
        )
    end
end

# ADDED: same limits in all four panels, then hide the repeated x labels
linkaxes!(vec(axs)...)                                 # CHANGED: vec() of the matrix
for r in 1:n_rows, c in 1:n_cols
    r < n_rows && hidexdecorations!(axs[r, c]; grid=false)   # x labels only on the bottom row
    c > 1      && hideydecorations!(axs[r, c]; grid=false)   # y labels only on the first column
end
colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

Colorbar(fig[1:n_rows, n_cols + 1], sc; label = "E [keV]")


Label(fig[0, 1:n_cols + 1],                                     # ADDED
    "Footpoints relative to $(station), ω/Ωₑ = $(ω_frac_used)";
    fontsize=16, font=:bold, tellwidth=false)

ylims!(axs[1], -0.8, -0.4)
#save("src/WPI/test_scripts/TOF/figures/footpoints_enu.png", fig)
