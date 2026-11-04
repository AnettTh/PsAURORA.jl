using AURORA
using CairoMakie
using Profile
using ProgressMeter
using JLD2


## Define magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Define the particles
μ = -cos(deg2rad(2))       # Almost field-aligned
L = 6.5
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))

E_grid = range(1e3, 40e3, length=60)
# Minimal energy of precipitation was 6.949 keV for this run

particle = [ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true) for E in E_grid]

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
plasma = PlasmaState(
    λ_grid,
    ϕ_eq,
    DentonDensity(magnetic_field, L, ϕ_eq),
    magnetic_field,
    Float64(L)
)

# Define the wave
Ω_e_eq = minimum(plasma.Ω_e)
ωs = range(Ω_e_eq*0.1, Ω_e_eq*0.5, length=60)

modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)   # 0.2 s chirp
amplitude = ConstantAmplitude(100e-12)                      # 100 pT
wave = WaveState(modulation, amplitude, plasma; Ψ0=0.0)

## Calculate TOF
TOF = zeros(length(E_grid), length(ωs))
use_boris = true

@showprogress for (i, p) in enumerate(particle)
    TOF[i, :] = WPI_TOF(
        ωs,
        p,
        plasma;
        wave=wave,
        use_boris=use_boris,
        field_dependent=true,
    )
end


outfile = "src/WPI/data/02112026TOF_heatmap.jld2"

jldsave(outfile;
    TOF,                              # [E, ω], size (60, 60)
    E_grid  = collect(E_grid),        # [eV]
    ωs      = collect(ωs),            # [rad/s]
    ω_frac  = collect(ωs) ./ Ω_e_eq,  # ω/Ωₑ,eq, the axis you plot against
    Ω_e_eq,
    μ, L, ϕ_eq,
    Bw      = amplitude.Bw,           # [T]
    Ψ0      = wave.Ψ0,
    t_chirp = modulation.t_duration,
    use_boris,
)

## Make Figure
data = load(outfile)

fig = Figure(size=(700, 500))
ax  = Axis(fig[1, 1];
    xlabel = "ω/Ωₑ",
    ylabel = "E [keV]",
    title  = "Time-of-flight using " * (data["use_boris"] ? "boris mover" : "integration") * "\n" *
             "α = $(round(rad2deg(acos(abs(data["μ"]))), digits=1))°" *
             ", MLT $(round(longitude_to_MLT(data["ϕ_eq"]; degrees=false), digits=1))" *
             ", L = $(data["L"])",
)

hm = heatmap!(ax, data["ω_frac"], data["E_grid"][10:end] ./ 1e3, data["TOF"][10:end, :]';
    colormap  = :turbo,
    nan_color = :white,
)
Colorbar(fig[1, 2], hm; label="TOF [s]")

#save("src/WPI/test_scripts/TOF/figures/02112026TOF_energy_heatmap.png", fig)
fig
