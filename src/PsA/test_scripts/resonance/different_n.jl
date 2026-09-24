using AURORA
using AURORA; c₀
using CairoMakie


# Shows that the ducted wave (θ=0°) can only scatter for the first and second harmonic
## Define particle and plasma
L = 6.5
E = 30e3
μ = - cos(deg2rad(1))
r0 = [L*RE, 0.0, 0.0]

particle = ParticleState(E, μ, r0, dipole_field; relativistic=true)

λ_grid = range(0.0, deg2rad(50), length=500)
#n_e0 = 1.8e7        # 18/cc, from Hsieh 2022
n_e0 = 5e6           # 5/cc, from Gan 2023 (random PhD) # TODO: figure out why this is too low to find anything??

plasma = PlasmaState(λ_grid, n_e0, dipole_field, L; ne_model=ne_denton)

# find and plot resonance
ns = [-1, 0, 1, 2]
ω_grid = range(plasma.Ω_e[1]*0.1, plasma.Ω_e[1]*0.5, length=500)
ω_grid_UB = range(plasma.Ω_e[1]*0.5, plasma.Ω_e[1]*0.9, length=500)

fig = Figure()
ax = Axis(
    fig[1,1];
    xlabel="ω [kHz]",
    ylabel="λ [°]",
    title="Resonance latitude as a function of wave frequency")

parameters_text = """
    n_e = $(n_e0/1e6) cm⁻³
    L = $L
    E = $(E/1e3) keV
    α = $(round(rad2deg(acos(abs(μ))), digits=1))°
    ω = 0.1 Ωₑ - 0.5 Ωₑ
    θ = 0°
    """

text!(ax, 0.02, 0.98;
    text       = parameters_text,
    space      = :relative,
    align      = (:left, :top),
    fontsize   = 11
)

for n in ns
    λ_res = resonance_latitude(ω_grid, particle, plasma; n=n)
    lines!(ax, ω_grid./1e3, rad2deg.(λ_res), label="n = $(Int(n))")
end

axislegend(ax, position=:rb)

##
save("src/PsA/test_scripts/resonance/figures/n_tests.png",fig)


##
fig = Figure(size=(1200, 500))

# Shared parameters text
parameters_text = """
    n_e = $(n_e0/1e6) cm⁻³
    L = $L
    E = $(E/1e3) keV
    α = $(round(rad2deg(acos(abs(μ))), digits=1))°
    θ = 0°
    """

# Panel 1: Lower band
ax1 = Axis(fig[1, 1];
    xlabel = "ω [kHz]",
    ylabel = "λ [°]",
    title  = "Lower-band chorus (0.1Ωₑ - 0.5Ωₑ)"
)

text!(ax1, 0.02, 0.98;
    text      = parameters_text * "    ω = 0.1Ωₑ - 0.5Ωₑ",
    space     = :relative,
    align     = (:left, :top),
    fontsize  = 11
)

for n in ns
    λ_res = resonance_latitude(ω_grid, particle, plasma; n=n)
    lines!(ax1, ω_grid ./ 1e3, rad2deg.(λ_res), label="n = $(Int(n))")
end

axislegend(ax1, position=:rb)

# Panel 2: Upper band
ax2 = Axis(fig[1, 2];
    xlabel = "ω [kHz]",
    ylabel = "λ [°]",
    title  = "Upper-band chorus (0.5Ωₑ - 0.9Ωₑ)"
)

text!(ax2, 0.02, 0.98;
    text      = parameters_text * "    ω = 0.5Ωₑ - 0.9Ωₑ",
    space     = :relative,
    align     = (:left, :top),
    fontsize  = 11
)

for n in ns
    λ_res = resonance_latitude(ω_grid_UB, particle, plasma; n=n)
    lines!(ax2, ω_grid_UB ./ 1e3, rad2deg.(λ_res), label="n = $(Int(n))")
end

axislegend(ax2, position=:rb)

# Link y-axes for easy comparison
linkyaxes!(ax1, ax2)

save("src/PsA/test_scripts/resonance/resonance_LB_UB.png", fig)
