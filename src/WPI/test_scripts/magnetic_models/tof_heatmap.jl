using AURORA
using CairoMakie
using Profile
using ProgressMeter
using JLD2


## Define magnetic field
ϕ = deg2rad(240.0)
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define the particles
μ = -cos(deg2rad(3))       # Almost field-aligned
L = 6.5
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

E_grid = range(10e3, 40e3, length=50)

particle = [ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true) for E in E_grid]

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
plasma = PlasmaState(λ_grid, ϕ_eq, DentonDensity(magnetic_field, L, ϕ_eq), magnetic_field, Float64(L))

# Define the wave
Ω_e_eq = minimum(plasma.Ω_e)
ωs = range(Ω_e_eq*0.1, Ω_e_eq*0.49, length=50)

# NOTE: I do not think this is correct anymore
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma)

## Calculate TOF
TOF = zeros(length(E_grid), length(ωs))
use_boris = true

@showprogress for (i, p) in enumerate(particle)
    TOF[i, :] = WPI_TOF(
        ωs,
        p,
        plasma,
        wave;
        use_boris=use_boris,
        field_dependent=true,
        #wave_launch_time=wave_chirp
    )
end


## Make Figure
fig = Figure(size=(700,500))
ax = Axis(
    fig[1,1],
    xlabel="ω/Ωₑ",
    ylabel="E [keV]",
    title="Time-of-flight using "
    * if use_boris "boris mover \n" else "integration \n " end
    * "α = $(round(rad2deg(acos(abs(μ))), digits=1))"
    * ", MLT $(round(longitude_to_MLT(ϕ; degrees=false), digits = 1))"
    * ", L=$(L)",
    #yscale=log10
)

hm = heatmap!(ax, ωs ./ Ω_e_eq, E_grid ./ 1e3, TOF; colormap=:turbo, nan_color=:white)
Colorbar(fig[1,2], hm; label="TOF [s]")
##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

save("src/WPI/test_scripts/magnetic_models/tof_heatmap_corrected.png", fig)
