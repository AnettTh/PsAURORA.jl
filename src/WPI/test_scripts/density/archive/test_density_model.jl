using AURORA
using CairoMakie
using Dates

##
L_vals = [4, 5, 6, 7, 8]
ϕ_eq = deg2rad(240.0)
time = DateTime("2020-01-01T02:00:00")
tsyganenko = TsyganenkoMagneticField(ϕ_eq; time=time)
dipole = DipoleMagneticField()
colors = [:blue, :red, :green, :orange, :purple]


##
fig = Figure(size=(1000, 400))

Label(fig[0, 1:2], "Dentons density model"; fontsize=18, font=:bold, tellwidth=false)

ax_tsyganenko = Axis(
    fig[1, 1];
    xlabel = "R [RE]",
    ylabel = "nₑ [m⁻³]",
    title  = "Adjusted model",
    yscale = log10
)
ax_dipole = Axis(
    fig[1, 2];
    xlabel = "R [RE]",
    ylabel = "nₑ [m⁻³]",
    title  = "Original model",
    yscale = log10
)

scatter!(ax_tsyganenko, [NaN], [NaN], label="(Lᵦ, Lₜₛ, Lₘₐₓ)", color=:transparent)
scatter!(ax_dipole, [NaN], [NaN], label="L", color=:transparent)

## Make Tsyganenko plot
for (L, color) in zip(L_vals, colors)
    model_tsyganenko = DentonDensity(tsyganenko, Float64(L), ϕ_eq; adjusted_model=true)
    no_adjustment_model_tsyganenko = DentonDensity(tsyganenko, Float64(L), ϕ_eq; adjusted_model=false)


    R_grid = range(RE + z_ionosphere, model_tsyganenko.R_max, length=500)

    n_e = [denton_density_model(
            Cartesian(R, 0.0, 0.0), model_tsyganenko.R_max, model_tsyganenko.L_Denton_bin
            ) * (1.0 + model_tsyganenko.duct_enhancement)
        for R in R_grid
    ]

    no_adjustment_n_e = [denton_density_model(
            Cartesian(R, 0.0, 0.0), no_adjustment_model_tsyganenko.R_max, no_adjustment_model_tsyganenko.L_Denton_bin
            ) * (1.0 + no_adjustment_model_tsyganenko.duct_enhancement)
        for R in R_grid
    ]


    lines!(
        ax_tsyganenko,
        collect(R_grid) ./ RE,
        n_e;
        color=color,
        label="($(round(model_tsyganenko.L_Denton_bin, digits=1)), " *
        "$(round(L, digits=1)), " *
        "$(round(model_tsyganenko.R_max / RE, digits=1)))"
    )
    lines!(
        ax_dipole,
        collect(R_grid) ./ RE,
        no_adjustment_n_e;
        color=color,
        label="$(round(L, digits=1))"
    )
end

##  Make dipole plot
for (L, color) in zip(L_vals, colors)
    model_dipole = DentonDensity(dipole, Float64(L), ϕ_eq)

    R_grid = range(RE + z_ionosphere, L * RE, length=500)

    n_e = [denton_density_model(
            Cartesian(R, 0.0, 0.0), model_dipole.R_max, model_dipole.L_Denton_bin
            ) * (1.0 + model_dipole.duct_enhancement)
        for R in R_grid
    ]

    lines!(
        ax_dipole,
        collect(R_grid) ./ RE,
        n_e;
        color=color,
        linestyle=:dash
    )
end
# Make labels to distinguish dipole and tsyganenko
lines!(ax_dipole, [NaN], [NaN], label="Tsyganenko", color=:black)
lines!(ax_dipole, [NaN], [NaN], label="Dipole", color=:black, linestyle=:dash)

## Set common scales for both models
linkyaxes!(ax_tsyganenko, ax_dipole)
hideydecorations!(ax_dipole; grid=false)
linkxaxes!(ax_tsyganenko, ax_dipole)
ticks = collect(1:0.5:9)
labels = [isinteger(t) ? string(Int(t)) : "" for t in ticks]
ax_tsyganenko.xticks = (ticks, labels)
ax_dipole.xticks = (ticks, labels)
axislegend(ax_tsyganenko, position=:lb)
axislegend(ax_dipole, position=:rt)
