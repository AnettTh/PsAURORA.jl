using AURORA
using CairoMakie


## Count how many of each energy and α₀ have precipitated
precipitated = @. isfinite(TOF) & (TOF > 0)
bounced = @. TOF == 0.0
resonant = @. !isnan(TOF)

collapse(A) = dropdims(sum(A; dims=(3, 4)); dims=(3, 4))

# Of size [E, α]
n_precip = collapse(precipitated)
n_bounced = collapse(bounced)
n_resonant = collapse(resonant)
n_total = size(TOF, 3) * size(TOF, 4)

probability = n_precip ./ n_total

## Make the figure
fig = Figure(size=(700, 300))
ax = Axis(
    fig[1,1],
    xlabel="α₀ [°]",
    ylabel="Probability",
    title="Scattering probability as a function of inital pitch angle \n $(params.:station) $(params.time)")

E_keV = collect(E_grid) ./ 1e3
α_deg = collect(α0_grid)

scatter!(ax, [NaN], [NaN], color=:transparent, label="Energy [keV]")
colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

for iE in eachindex(E_keV)
    scatter!(ax, α_deg, probability[iE, :]; color=colors[iE], label="$(E_keV[iE])")
end

scatter!(ax, [NaN], [NaN], color=:transparent, label="$(length(TOF)) total")
scatter!(ax, [NaN], [NaN], color=:transparent, label="$(sum(n_precip)) scattered")
scatter!(ax, [NaN], [NaN], color=:transparent, label="$(sum(n_bounced)) bounced")
scatter!(ax, [NaN], [NaN], color=:transparent, label="$(sum(n_resonant) - length(TOF)) non-resonant")

Legend(fig[1, 2], ax)
xlims!(ax, 0, 10.1)
save("src/WPI/test_scripts/big_data/figures/11202026_ScatteringProbability.png", fig)
