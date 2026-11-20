using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Configuration
duct_enhancement = 0.0
station = :sodankyla
n_T = 100_000

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
E = 30e3
α0 = 0.5
μ0 = -cos(deg2rad(α0))
particle = ParticleState(E, μ0, fieldline)

## Define the density and plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)

Ψ0_grid = range(0, 2π, length=5)[1:end-1]

waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Frequencies to evaluate
ω_frac = 0.3
ω =  ω_frac * minimum(plasma.Ω_e)

## Run Boris-mover
results = Vector{Any}(nothing, length(Ψ0_grid))

@showprogress for i in eachindex(waves)
    # TODO: Add r_resonance to run boris_mover
    #results[i] = boris_mover(particle, fieldline, waves[i]; store_trajectory=true, stop_when_bouncing=false, n_T=n_T, ω=ω)
    results[i] = boris_mover(
        particle;
        wave = waves[i],
        store_trajectory = true,
        stop_when_bouncing = false,
        n_T = n_T
    )
end

## Latitude and pitch angle of each trajectory
function λ_and_α(r)
    n = findlast(>(0), r.time)
    x, y, z = r.position[1:n, 1], r.position[1:n, 2], r.position[1:n, 3]
    λ = rad2deg.(asin.(z ./ sqrt.(x.^2 .+ y.^2 .+ z.^2)))
    return λ, rad2deg.(r.α[1:n])
end

# Adiabatic pitch angle along the line, without a wave (does not depend on Ψ₀)
α_ad = [begin
            a = pitch_angle(particle.α_eq, particle.B_eq, B)
            isnothing(a) ? NaN : rad2deg(a)
        end for B in fieldline.Bs]
λ_line = rad2deg.(fieldline.λs)

## Figure: one axis per Ψ₀
n_rows, n_cols = 2, 2                                  # needs n_rows * n_cols == length(Ψ0_grid)
fig = Figure(size=(800, 400))
panel_axes = Matrix{Axis}(undef, n_rows, n_cols)

λ_start, _ = λ_and_α(results[1])                       # same start point for every phase

for (i, r) in enumerate(results)
    row, col = fldmod1(i, n_cols)
    ax = Axis(fig[row, col];
        xlabel = "λ (SM) [°]",
        ylabel = "α [°]",
    )
    panel_axes[row, col] = ax

    λ, α = λ_and_α(r)
    lines!(ax, λ, α; color=colors[i])
    lines!(ax, λ_line, α_ad; color=:black, linestyle=:dash)
    vlines!(ax, [λ_start[1]]; color=:black, linestyle=:dot)
    if !isnothing(r.footpoint)
        scatter!(ax, [λ[end]], [α[end]];
            color=colors[i], marker=:star5, markersize=18,
            strokecolor=:black, strokewidth=1)
        println("precipitated")
    end
end

ax_leg = panel_axes[1, 1]
[lines!(ax_leg, [NaN], [NaN]; color=colors[i], label="Ψ₀ = $(round(Ψ0_grid[i] / π, digits=1)) π") for i in eachindex(results)]
lines!(ax_leg, [NaN], [NaN]; color=:black, linestyle=:dash, label="Adiabatic")
lines!(ax_leg, [NaN], [NaN]; color=:black, linestyle=:dot, label="λᵣₑₛ")
lines!(ax_leg, [NaN], [NaN]; color=:transparent, label="E = $(E / 1e3) keV")
lines!(ax_leg, [NaN], [NaN]; color=:transparent, label="α₀ = $(α0)°")
lines!(ax_leg, [NaN], [NaN]; color=:transparent, label="n_T = $(n_T)")

linkaxes!(vec(panel_axes)...)
for row in 1:n_rows, col in 1:n_cols
    row < n_rows && hidexdecorations!(panel_axes[row, col]; grid=false)
    col > 1      && hideydecorations!(panel_axes[row, col]; grid=false)
end
Legend(fig[1:n_rows, n_cols + 1], ax_leg;
    framevisible = false,
    halign       = :left,
    valign       = :top,
    labelsize    = 12,
)

Label(fig[0, 1:n_cols], "Pitch angle evolution along particle trajectory \n$station $t";
    fontsize=16, font=:bold, tellwidth=false)

colgap!(fig.layout, 30)
rowgap!(fig.layout, 20)

##
#save("src/WPI/test_scripts/pitch_angle/figures/11202026_AlphaLatitudeEvolutionGrid.png", fig)
#
#ylims!(panel_axes[1, 1], 0, 20)
#xlims!(panel_axes[1, 1], -5, 25)
#
#save("src/WPI/test_scripts/pitch_angle/figures/11202026_AlphaLatitudeEvolutionGrid_zoomed.png", fig)
