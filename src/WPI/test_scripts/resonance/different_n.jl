using AURORA
using AURORA; c₀
using CairoMakie


# Shows that the ducted wave (θ=0°) can only scatter for the first and second harmonic
## Define particle and plasma
L = 6.5
E = 30e3
μ = - cos(deg2rad(3))
r0 = [L*RE, 0.0, 0.0]

particle = ParticleState(E, μ, r0, dipole_field)

λ_grid = range(0.0, deg2rad(50), length=500)
#n_e0 = 5e6        # 18/cc, from Hsieh 2022

plasma = PlasmaState(λ_grid, ne_denton, dipole_field, L)

# find and plot resonance
ns = [-1, 0, 1, 2]
ω_grid = range(plasma.Ω_e[1]*0.1, plasma.Ω_e[1]*0.5, length=500)

##
fig = Figure()
ax = Axis(
    fig[1,1];
    xlabel="ω [kHz]",
    ylabel="λ [°]",
    title="Resonance latitude as a function of wave frequency")

parameters_text = """
    n_e = Denton-model
    L = $L
    E = $(E/1e3) keV
    α = $(round(rad2deg(acos(abs(μ))), digits=1))°
    ω = 0.25 Ωₑ - 0.5 Ωₑ
    θ = 0°
    """

text!(ax, 0.02, 0.98;
    text       = parameters_text,
    space      = :relative,
    align      = (:left, :top),
    fontsize   = 11
)
##
for n in ns
    λ_res = resonance_latitude(ω_grid, particle, plasma; n=n)
    lines!(ax, ω_grid./1e3, rad2deg.(λ_res), label="n = $(Int(n))")
end

axislegend(ax, position=:rb)

##
save("src/WPI/test_scripts/resonance/n_tests.png",fig)
