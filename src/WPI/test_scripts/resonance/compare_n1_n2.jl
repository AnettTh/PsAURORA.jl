using AURORA
using AURORA; c₀
using CairoMakie
using ProgressMeter

## Make magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define the particles
μ = -cos(deg2rad(2))       # Almost field-aligned
L = 6.5
r0 = Cartesian(L*RE, 0.0, ϕ_eq)

E_grid = range(10e3, 40e3, length=100)

particles = @showprogress [ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true) for E in E_grid]

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=500)
plasma = PlasmaState(λ_grid, DentonDensity(magnetic_field, L, ϕ_eq))


## Define the wave
Ωe_frac = range(0.1, 0.5, length=5)
ωs = plasma.Ω_e[1] .* Ωe_frac

modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

##
TOF_n1 = zeros(length(E_grid), length(ωs))
TOF_n2 = zeros(length(E_grid), length(ωs))
## NOTE: This does not run anymore, as i probably have fixed some bug and n=2 now is out of valid range
@showprogress for (i, p) in enumerate(particles)
    TOF_n1[i, :] = WPI_TOF(ωs, p, plasma; wave, n=1, use_boris=true)
    TOF_n2[i, :] = WPI_TOF(ωs, p, plasma; wave, n=2, use_boris=true)
end

# NOTE: The resonance here is at positive latitudes. Is this some old hardcoded bug or due to n=2 being used?

##
fig = Figure(size=(1200, 500))

Label(fig[0, 1:3],
    "Wave-particle interaction time-of-flight";
    fontsize=16, font=:bold
)

ax1 = Axis(fig[1, 1];
    xlabel = "Time [s]",
    ylabel = "E [keV]",
    title  = "n = 1"
)

ax2 = Axis(fig[1, 2];
    xlabel = "Time [s]",
    ylabel = "E [keV]",
    title  = "n = 2"
)

colors = [:blue, :red, :green, :orange, :purple]

for (i, ω) in enumerate(ωs)
    label = "ω = $(round(ω/plasma.Ω_e[1], digits=2))Ωe"
    lines!(ax1, TOF_n1[:, i], E_grid ./ 1e3; color=colors[i], label=label)
    lines!(ax2, TOF_n2[:, i], E_grid ./ 1e3; color=colors[i], label=label)
end

params = """
    Parameters:
    L = $L
    α = $(round(rad2deg(acos(abs(μ))), digits=1))°
    n_e = denton-model
    """

Label(fig[1, 3],
    params;
    fontsize   = 10,
    tellwidth  = false,
    tellheight = false,
    halign     = :left,
    valign     = :bottom
)

# Link axes so they share the same limits
linkyaxes!(ax1, ax2)
linkxaxes!(ax1, ax2)

Legend(fig[1, 3], ax1)
