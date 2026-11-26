using AURORA
using JLD2
using CairoMakie
using Statistics
using LinearAlgebra

## Load and select the pitch angle
data = load("src/WPI/data/11-25-2026/energy-time-position.jld2")

## Get the data and reshape it
iα = 3
α = data["α_grid"][iα]                    # [°]
TOF_all = data["TOF"][:, iα, :, :]        # [E, μ, ω, Ψ₀]
dx_all = data["fp_enu"][1, :, iα, :, :]   # [xyz, E, μ, ω, Ψ₀]
E_keV = data["E_grid"] ./ 1e3             # [keV]

E_all = repeat(E_keV, 1, size(dx_all, 2), size(dx_all, 3))

# Pick out only valid cases
valid = @. isfinite(dx_all) & isfinite(TOF_all) & (TOF_all > 0)

# Flat vectors using the same masks
t = TOF_all[valid] .- minimum(TOF_all[valid])
dx = dx_all[valid] ./ 1e3
E = E_all[valid]

## Do statistics on the dataset relative to the energy
function quantile_per_energy(values, E_values, E_grid, q)
    # Resutl is NaN where there is no valid samples in the energy
    result = fill(NaN, length(E_grid))

    # Loop over all energies
    for (iE, E_i) in enumerate(E_grid)
        samples = values[E_values .== E_i]
        isempty(samples) && continue
        result[iE] = quantile(samples, q)
    end
    return result
end

dx_med = quantile_per_energy(dx, E, E_keV, 0.5)
t_min  = quantile_per_energy(t , E, E_keV, 0.0)
t_q01  = quantile_per_energy(t , E, E_keV, 0.01)
t_q05  = quantile_per_energy(t , E, E_keV, 0.05)

function fit_poly(x, y, deg)
    # Get limits
    ok = isfinite.(y)
    x0, x1 = extrema(x[ok])

    # Scale to [-1, 1], keeping the fit well-conditioned
    s(xi) = 2 * (xi - x0) / (x1 - x0) - 1
    X = [s(xi)^k for xi in x[ok], k in 0:deg]
    c = X \ y[ok]
    return xi -> sum(c[k+1] * s(xi)^k for k in 0:deg)
end

## Make polynomial fit
dx_fit = fit_poly(E_keV, dx_med, 3)
tmin_fit = fit_poly(E_keV, t_min, 3)
tq01_fit = fit_poly(E_keV, t_q01, 3)
tq05_fit = fit_poly(E_keV, t_q05, 3)

##
τmin = t .- tmin_fit.(E)
τq01 = t .- tq01_fit.(E)
τq05 = t .- tq05_fit.(E)

E_groups = [(7, 12), (12, 20), (20, 30), (30, 40)]    # keV
group_colors = ["#E69F00", "#56B4E9", "#009E73", "#D55E00"]

fig = Figure(size=(900, 550))
#Label(fig[0, 1], text="Arrival time compared to onset curve", font=:bold, fontsize=12)

ax1 = Axis(fig[1, 1]; ylabel="# Events", title="Arrival time compared to onset curve, tₘᵢₙ")
for ((lo, hi), c) in zip(E_groups, group_colors)
    m = (E .≥ lo) .& (E .< hi)
    count(m) < 5 && continue
    stephist!(ax1, τmin[m]; bins=30, normalization=:pdf, color=c, linewidth=2, label="$lo-$hi keV")
end

ax2 = Axis(fig[2, 1]; xlabel="τ = t − t₀(E) [s]", ylabel="# Events", title="t₀₁")

for ((lo, hi), c) in zip(E_groups, group_colors)
    m = (E .≥ lo) .& (E .< hi)
    count(m) < 5 && continue
    stephist!(ax2, τq01[m]; bins=30, normalization=:pdf, color=c, linewidth=2, label="$lo-$hi keV")
end

ax3 = Axis(fig[3, 1]; xlabel="τ = t − t₀(E) [s]", ylabel="# Events", title="t₀₅")

for ((lo, hi), c) in zip(E_groups, group_colors)
    m = (E .≥ lo) .& (E .< hi)
    count(m) < 5 && continue
    stephist!(ax3, τq05[m]; bins=30, normalization=:pdf, color=c, linewidth=2, label="$lo-$hi keV")
end

axislegend(ax1)

linkxaxes!(ax1, ax2, ax3)
hidexdecorations!(ax1, grid = false)
hidexdecorations!(ax2, grid = false)

save("src/WPI/test_scripts/input_flux/figures/11262026_OnsetCurveHistogram.png", fig)
