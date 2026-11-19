using AURORA
using CairoMakie
using ProgressMeter

##
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T03:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    t;
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)

ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define the particle
α_deg = 0.5
μ0 = -cos(deg2rad(α_deg))
E_grid = range(7e3, 40e3, length=100)

particles = [ParticleState(E, μ0, fieldline) for E in E_grid]

## Define the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(density, fieldline)

## Define the wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)
wave = WaveState(modulation, amplitude, plasma, fieldline)

ω_frac = range(0.1, 0.5, length=5)
ωs = ω_frac * minimum(plasma.Ω_e)

## Make figure
colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

fig = Figure(size=(600, 400))
ax  = Axis(fig[1, 1];
    xlabel = "E [keV]",
    ylabel = "λᵣₑₛ [°]",
    title  = "Resonance latitude vs energy \n" *
    if sign(μ0) < 0 "propagating opposite to the wave" else "propagating in the same direction as the wave" end,
    #xscale = log10
)

##
for (ω, color) in zip(ω_frac, colors)

    λ_res_grid = @showprogress map(E_grid) do E
        p = ParticleState(E, μ0, r0, magnetic_field, ϕ_eq; relativistic=true)
        λ_res = point_of_resonance(ω, p, plasma; n=n).λ
        isnan(λ_res) ? NaN : rad2deg(λ_res)
    end

    lines!(ax, collect(E_grid)./1e3, λ_res_grid; color=color, label="ω = $ω Ωₑ")
end

scatter(ax, [NaN], [NaN], color=:transparent, label="α₀ = $(α0)°")
scatter(ax, [NaN], [NaN], color=:transparent, label="n = $n°")

axislegend(ax, position=:rb)

#save("src/WPI/test_scripts/resonance/figures/11192026_ResonanceLatitudeEnergy.png", fig)
