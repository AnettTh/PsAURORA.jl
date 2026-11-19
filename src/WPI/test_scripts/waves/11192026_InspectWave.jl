using AURORA
using CairoMakie
using Dates

##
duct_enhancement = 0.0
station = :sodankyla

## Field, station line, particle and plasma
t = DateTime("2020-12-21T02:00:00")
magnetic_field = TsyganenkoMagneticField(t; pdyn=2.0, dst=-28.0, byimf=0.0, bzimf=2.0)
fieldline = StationFieldLine(magnetic_field, LTPReference(station), t)

E = 30e3
α0 = 0.5
μ0 = -cosd(α0)
particle = ParticleState(E, μ0, fieldline)

density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
amplitude = ConstantAmplitude(100e-12)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WaveState(modulation, amplitude, plasma, fieldline)


## Make a wave field
time = 0.0:0.001:0.2     # Current time
s = 1.0                  # Current position
Ψ = 0.0                  # Current wave phase

B0 = particle.B0    # Initial position magnetic field
λ        = fieldline.λ_eq          # use the latitude of B0, not 0.0
s        = 0.0                     # only matters for position-dependent amplitudes
_, e1, e2, _ = magnetic_basis(B0)

## Time axis
ts = collect(0.0:0.00001:0.2)        # [s]

## Phase seen by a particle at rest: Ψ(t) = Ψ0 - ∫ ω dt
ω_t = [wave_frequency(wave, τ) for τ in ts]
Ψ_t = wave.Ψ0 .- cumsum([0.0; diff(ts) .* ω_t[2:end]])

## Evaluate the fields
out = [wave_fields(wave, τ, s, Ψ, B0, λ) for (τ, Ψ) in zip(ts, Ψ_t)]

δB1 = [dot(o.δB, e1) for o in out] .* 1e12      # [pT]
δB2 = [dot(o.δB, e2) for o in out] .* 1e12
δE1 = [dot(o.δE, e1) for o in out] .* 1e3       # [mV/m]
δE2 = [dot(o.δE, e2) for o in out] .* 1e3
k_par = [o.k_parallel for o in out]
Ω_eq  = minimum(plasma.Ω_e)

## Figure
fig = Figure(size=(900, 400))

ax1 = Axis(fig[1, 1]; xlabel="t [s]", ylabel="ω / Ωₑ", title="Wave frequency")
lines!(ax1, ts, ω_t ./ Ω_eq)

ax2 = Axis(fig[1, 2]; xlabel="t [s]", ylabel="k∥ [rad/km]", title="Wave number at the evaluation point")
lines!(ax2, ts, k_par * 1e3)

ax3 = Axis(fig[2, 1]; xlabel="t [s]", ylabel="δB [pT]", title="δB components along e1, e2")
lines!(ax3, ts, δB1; label="δB·e1")
lines!(ax3, ts, δB2; label="δB·e2")
axislegend(ax3)

ax4 = Axis(fig[2, 2]; xlabel="t [s]", ylabel="δE [mV/m]", title="δE components along e1, e2")
lines!(ax4, ts, δE1; label="δE·e1")
lines!(ax4, ts, δE2; label="δE·e2")
axislegend(ax4)

rowsize!(fig.layout, 1, Relative(1/3))
rowsize!(fig.layout, 2, Relative(2/3))

linkxaxes!(ax1, ax3)
linkxaxes!(ax2, ax4)
hidexdecorations!(ax1; grid=false)
hidexdecorations!(ax2; grid=false)

save("src/WPI/test_scripts/waves/figures/11192026_InspectWave.png")
