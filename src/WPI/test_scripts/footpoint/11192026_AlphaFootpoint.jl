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
E = range(7e3, 40e3, length=1000)
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
ω = ω_frac * minimum(plasma.Ω_e)

α_footpoint = fill(NaN, length(E))
FP = fill(NaN, 3, length(E))
progress = Progress(length(particle))
Threads.@threads :greedy for i in eachindex(particle)
    result = boris_mover(particle[i], fieldline, wave; ω=ω)
    if !isnothing(result.footpoint)
        α_footpoint[i] = result.α_footpoint
        FP[:, i] .= result.footpoint
    end
    next!(progress)
end

fp_enu = fill(NaN, size(FP))
for iE in axes(FP, 2)
    any(isnan, @view FP[:, iE]) && continue
    fp_enu[:, iE] .= FP_to_ENU(FP[:, iE], LTPReference(station), t)
end

east_km = fp_enu[1, :] ./ 1e3

E_keV = collect(E) ./ 1e3
##
fig = Figure(size=(700, 400))
ax = Axis(
    fig[1, 1],
    xlabel="E [keV]",
    ylabel="α (footpoint) [°]",
    title="Pitch angle at footpoint for precipitating particles at $station $t \n" * "α₀ = $α0"
)

lim = maximum(abs, filter(isfinite, east_km); init=1.0)
α_fp_plot = [isnothing(a) ? NaN : rad2deg(a) for a in α_footpoint]
sc = scatter!(ax, E_keV, α_fp_plot, color=east_km, colormap=:turbo, colorrange=(0, lim))
Colorbar(fig[1, 2], sc; label="East drift [km]")
ylims!(ax, 0, 90)
xlims!(ax, 7, 30)

save("src/WPI/test_scripts/footpoint/figures/11192026_AlphaFootpoint.png", fig)
