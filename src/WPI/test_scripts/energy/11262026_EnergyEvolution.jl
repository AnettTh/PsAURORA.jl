using AURORA
using CairoMakie
using ProgressMeter
using Dates

## Define configuration
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
ref = LTPReference(station)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Define the particles
α0 = 1.5
μ0 = -cos(deg2rad(α0))
E_list = [7e3, 15e3, 30e3, 40e3]
particles = [ParticleState(E, μ0, fieldline) for E in E_list]

## Define the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define a wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)
Ψ0_grid = range(0.0, 2π, length=7)[1:end-1]

waves = [WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0) for Ψ0 in Ψ0_grid]

## Run Boris-mover (nothing = no resonance for that particle)
results = Matrix{Any}(nothing, length(particles), length(Ψ0_grid))   # [energy, phase]

@showprogress for iE in eachindex(particles), i in eachindex(Ψ0_grid)
    results[iE, i] = boris_mover(particles[iE], fieldline, waves[i]; store_trajectory=true)
end

## Make colorscale
colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

## Energy-change plot, one axis per energy
n_rows, n_cols = 2, 2
fig = Figure(size=(700, 450))
Label(fig[0, 1:2], text="Change in energy along particle trajectory for α=$(α0)°", font=:bold, fontsize=16, tellwidth = false)
panel_axes = Matrix{Axis}(undef, n_rows, n_cols)

for k in eachindex(E_list)
    row, col = fldmod1(k, n_cols)
    ax = Axis(fig[row, col];
        xlabel = "t [s]",
        ylabel = "ΔE [%]",
    )
    panel_axes[row, col] = ax
    energy_color = colors[mod1(k, length(colors))]

    for i in eachindex(Ψ0_grid)
        isnothing(results[k, i]) && continue
        E0 = results[k, i].E[1]
        lines!(ax, results[k, i].time, ((results[k, i].E .- E0) ./ E0) .* 100;
            color=energy_color)
    end
    lines!(ax, [NaN], [NaN], color=energy_color, label="$(E_list[k] / 1e3) keV")
    axislegend(ax)
end

## Link x within each column and y within each row
for col in 1:n_cols
    linkxaxes!(panel_axes[:, col]...)
end
for row in 1:n_rows
    linkyaxes!(panel_axes[row, :]...)
end

## Hide repeated ticks and labels, keep the grid lines
for row in 1:n_rows, col in 1:n_cols
    row < n_rows && hidexdecorations!(panel_axes[row, col]; grid=false)
    col > 1 && hideydecorations!(panel_axes[row, col]; grid=false)
end

xlims!(panel_axes[1, 1], 0.0, 1.5)
xlims!(panel_axes[1, 2], 0.0, 1.5)
ylims!(panel_axes[1, 1], -0.2, 0.8)
ylims!(panel_axes[2, 1], -0.2, 0.8)

#legend_elements = [LineElement(color=colors[mod1(k, length(colors))], linewidth=2)
#                   for k in eachindex(E_list)]
#legend_labels = ["$(E / 1e3) keV" for E in E_list]
#
#axislegend(panel_axes[2, 2])
#Legend(fig[1:n_rows, n_cols + 1], legend_elements, legend_labels, "Energy"; framevisible=false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 30)

save("src/WPI/test_scripts/energy/figures/11262026_EnergyEvolution.png", fig)
