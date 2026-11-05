using AURORA
using CairoMakie
using ProgressMeter
using JLD2
using Dates


## Define method
name = "src/WPI/data/11-03-2026/TOF_duct.jld2"
use_boris = true
relativistic = true
duct_enhancement = 0.5
WPI = true
field_dependent = true
n = 1

## Define magnetic field
ϕ_eq = deg2rad(240.0)
t = DateTime("2020-01-01T02:00:00")
pdyn = 2.0
dst = -50.0
byimf = 1.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    ϕ_eq;
    time=t,
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)

## Define the particles
μ = -cos(deg2rad(0.5))
L = 6.5
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))
E_grid = range(7e3, 40e3, length=20)

particles = Vector{ParticleState}(undef, length(E_grid))
@showprogress for (i, E) in enumerate(E_grid)
    particles[i] = ParticleState(
        E,
        μ,
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=relativistic,
    )
end

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=1000)
lb_low, lb_high = 0.1, 0.5
density = DentonDensity(magnetic_field, L, ϕ_eq; duct_enhancement=duct_enhancement)

plasma = PlasmaState(λ_grid, density; lb_low=lb_low, lb_high=lb_high)

## Define the wave
ωs = range(plasma.ω_lb[1], plasma.ω_lb[end], length=5)
ω_frac = range(lb_low, lb_high, length=length(ωs))
amplitude = ConstantAmplitude(100e-12)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WPI ? WaveState(modulation, amplitude, plasma; wave_direction=-1) : nothing


## Calculate TOF
TOF = zeros(length(E_grid), length(ωs))

progress = Progress(length(particles))
Threads.@threads for i in eachindex(particles)
    TOF[i, :] = WPI_TOF(
        ωs,
        particles[i],
        plasma;
        wave=wave,
        use_boris=use_boris,
        field_dependent=field_dependent,
        n=n
    )
    next!(progress)
end


## Save the data
magnetic_field_model = magnetic_field isa TsyganenkoMagneticField ? "Tsyganenko" : "Dipole"
parameters = (
    L = L,
    μ = μ,
    ϕ_eq = ϕ_eq,
    ω_frac = ω_frac,
    magnetic_field_model = magnetic_field_model,
    t = t,
    pdyn = pdyn,
    dst = dst,
    byimf = byimf,
    bzimf = bzimf,
    ps = ps,
    use_boris = use_boris,
    relativistic = relativistic,
    duct_enhancement = duct_enhancement,
    WPI = WPI,
    field_dependent = field_dependent,
    n = n
)

@save name TOF E_grid ωs parameters



## Make the figure
fig = Figure(size=(900, 500))
ax = Axis(
    fig[1,1],
    xlabel="TOF [s]",
    ylabel="E [keV]",
    title="TOF at $(longitude_to_MLT(ϕ_eq; degrees=false)) MLT with α = " *
    "$(round(rad2deg(acos(abs(μ))), digits=1))°",
)

colors = [:blue, :red, :green, :orange, :purple]  # one per ω


for (i, ω) in enumerate(ωs)
    tof = replace(TOF[:, i], 0.0 => NaN)
    scatter!(ax, tof, E_grid ./ 1e3; color = colors[i], label="ω/Ωₑ=$(ω_frac[i])")
end

information = [
    "L = $(L)",
    "$(magnetic_field_model) field model",
    "$(t)",
    "pdyn = $(pdyn)",
    "dst = $(dst)",
    "byimf = $(byimf)",
    "bzimf = $(bzimf)",
    "ps = $(ps)",
    "use_boris = $(use_boris)",
    "relativistic = $(relativistic)",
    "$(duct_enhancement) duct enhancement",
    "WPI = $(WPI)",
    "field_dependent = $(field_dependent)",
    "n = $(n)",
]
for i in eachindex(information)
    lines!(ax, [NaN], [NaN], label=information[i], color=:transparent)
end

Legend(fig[1,2], ax)

##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

#save("src/WPI/test_scripts/energies/tsyg_tof.png", fig)
#save("src/WPI/test_scripts/TOF/figures/TOF.png", fig)
