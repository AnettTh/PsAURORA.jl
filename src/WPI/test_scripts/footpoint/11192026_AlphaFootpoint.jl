using AURORA
using Dates
using ProgressMeter
using CairoMakie

## Configuration
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define particle
E = range(7e3, 30e3, length=1000)
α0 = 1.0
μ0 = -cos(deg2rad(α0))

particle = [ParticleState(e, μ0, fieldline) for e in E]

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

wave = WaveState(modulation, amplitude, plasma, fieldline)

## Frequencies to evaluate
ω_frac = 0.3
ω =  ω_frac * minimum(plasma.Ω_e)

α_footpoint = fill(NaN, length(E))
@showprogress for (i, p) in enumerate(particle)
    result = boris_mover(p; wave=wave)
    isnothing(result.footpoint) && continue
    α_footpoint[i] = result.α_footpoint
end

E = collect(E) ./ 1e3

fig = Figure(size=(700, 400))
ax = Axis(
    fig[1, 1],
    xlabel="E [keV]",
    ylabel="α (footpoint) [°]",
    title="Pitch angle at footpoint for precipitating particles at $station $t \n" * "α₀ = $α0"
)

α_fp_plot = [isnothing(a) ? NaN : rad2deg(a) for a in α_footpoint]
sc = scatter!(ax, E, α_fp_plot, color="#E69F00")
ylims!(ax, 0, 90)
xlims!(ax, 7, 30)

save("src/WPI/test_scripts/footpoint/figures/11192026_AlphaFootpoint.png", fig)
