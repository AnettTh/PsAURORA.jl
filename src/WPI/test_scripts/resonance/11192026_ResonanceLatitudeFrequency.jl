using AURORA
using CairoMakie
using Dates

## Define model
ns = [1, 2]
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T03:00:00")
ref = LTPReference(:sodankyla)
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define the particle
α_deg = 0.5
μ0 = -cos(deg2rad(α_deg))
E_grid = [7e3, 15e3, 25e3, 40e3]

particles = [ParticleState(E, μ0, fieldline) for E in E_grid]

## Define the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma, fieldline)

ω_grid = range(0.1, 0.5, length=100) * minimum(plasma.Ω_e)

## Do computation
results = Dict{Tuple{Int, Int}, Any}()
for n in ns, iE in eachindex(E_grid)
    results[(n, iE)] = point_of_resonance(ω_grid, particles[iE], fieldline, plasma, wave; n=n)
end

## Make figure
fig = Figure(size = (800, 400))
ax = Axis(
    fig[1,1];
    xlabel="ω / Ωₑ",
    ylabel="λ [°]",
    title="Resonance latitude as a function of wave frequency")

colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
linestyles = [:solid, :dash]
Ω_eq = minimum(plasma.Ω_e)

## Plot resonances
for n in ns, iE in eachindex(E_grid)
    res = results[(n, iE)]
    λ_res  = [isnothing(r) ? NaN : r.λ for r in res]
    Vr_res = [isnothing(r) ? NaN : r.v for r in res]
    lines!(
        ax,
        ω_grid / Ω_eq,
        rad2deg.(λ_res);
        color = colors[iE],
        linestyle = linestyles[n],
    )
end

legend_elements = vcat(
    [LineElement(color = colors[iE]) for iE in eachindex(E_grid)],
    [LineElement(color = :black, linestyle = linestyles[n]) for n in ns],
)
legend_labels = vcat(
    ["E = $(round(E / 1e3, digits=1)) keV" for E in E_grid],
    ["n = $(n)" for n in ns],
)

## Add sub-grid
side = GridLayout(fig[1, 2]; valign = :top)

##
config_info = join([
    "Station: $(station)",
    "Coorinates: $(LTPReference(station).lat_deg)°, $(LTPReference(station).lon_deg)°",
    "Time: $(t)",
    "R_eq = $(round(fieldline.r_eq / RE, digits=2))",
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
    side[2, 1],
    legend_elements,
    legend_labels;
    framevisible = false,
    halign       = :left,
    labelsize    = 11,
)

rowgap!(side, 12)

##
save("src/WPI/test_scripts/resonance/figures/11192026_ResonanceLatitudeFrequency.png", fig)
