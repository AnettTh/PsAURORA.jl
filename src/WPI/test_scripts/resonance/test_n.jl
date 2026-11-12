using AURORA
using CairoMakie


## Define model
use_boris        = true
field_dependent  = true
relativistic     = true
WPI              = true
ns               = [1, 2]
duct_enhancement = 0.0
station          = :sodankyla

## Define magnetic field
t     = DateTime("2020-12-21T03:00:00")
ref   = LTPReference(:sodankyla)
pdyn  = 2.0
dst   = -28.0
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
α_deg = 0.5
μ0 = -cos(deg2rad(α_deg))

trace_L_shell = find_L_shell(magnetic_field, ref, t)
L = trace_L_shell.L
position_eq = trace_L_shell.position_eq

ϕ_eq = magnetic_field.ϕ
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))
Es = [7e3, 15e3, 25e3, 40e3]

particles = [
    ParticleState(E, μ0, r0, magnetic_field, ϕ_eq; relativistic=relativistic)
    for E in Es
]

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
lb_low, lb_high = 0.1, 0.5
density = DentonDensity(
    magnetic_field,
    L,
    ϕ_eq;
    duct_enhancement=duct_enhancement,
)

plasma = PlasmaState(λ_grid, density; lb_low=lb_low, lb_high=lb_high)

## Define the wave
ω_grid = plasma.ω_lb

## Do computation
results = Dict{Tuple{Int, Int}, Any}()
for n in ns, iE in eachindex(Es)
    results[(n, iE)] = point_of_resonance(ω_grid, particles[iE], plasma; n=n)
end

## Make figure
fig = Figure(size = (800, 400))
ax = Axis(
    fig[1,1];
    xlabel="ω / Ωₑ",
    ylabel="λ [°]",
    title="Resonance latitude as a function of wave frequency")

colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
linestyles = Dict(1 => :solid, 2 => :dash)
Ω_eq = minimum(plasma.Ω_e)      # TODO: add this to plasmastate!

## Plot resonances
for n in ns, iE in eachindex(Es)
    λ_res, Vr_res = results[(n, iE)]
    lines!(
        ax,
        ω_grid / Ω_eq,
        rad2deg.(λ_res);
        color     = colors[iE],
        linestyle = linestyles[n],
    )
end

legend_elements = vcat(
    [LineElement(color = colors[iE]) for iE in eachindex(Es)],
    [LineElement(color = :black, linestyle = linestyles[n]) for n in ns],
)
legend_labels = vcat(
    ["E = $(round(E / 1e3, digits=1)) keV" for E in Es],
    ["n = $(n)" for n in ns],
)

## Add sub-grid
side = GridLayout(fig[1, 2]; valign = :top)

##
config_info = join([
    "Station: $(station)",
    "Coorinates: $(LTPReference(station).lat_deg)°, $(LTPReference(station).lon_deg)°",
    "Time: $(t)",
    "L = $(round(L, digits=2))",
    "",
    "pdyn = $(pdyn) nPa",
    "Dst = $(dst) nT",
    "By = $(byimf) nT",
    "Bz = $(bzimf) nT",
    "tilt = $(round(rad2deg(magnetic_field.ps), digits=1))°",
    "",
    "Pitch-angle: $(α_deg)°"
], "\n")

Label(
    side[1, 1],
    config_info;
    justification = :left,        # left-aligned lines inside the box
    halign        = :left,
    valign        = :top,
    tellheight    = true,        # the box does not set the height of the rows
    fontsize      = 11,
)
Legend(
    side[2, 1],                     # below the info text
    legend_elements,
    legend_labels;
    framevisible = false,           # no box
    halign       = :left,
    labelsize    = 11,
)

rowgap!(side, 12)

##
save("src/WPI/test_scripts/resonance/figures/tsyg_n_tests.png",fig)
