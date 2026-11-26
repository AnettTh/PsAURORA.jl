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

## Plots
fig = Figure(size=(750, 350))

ax1 = Axis(fig[1, 1]; xlabel="dx [km]", ylabel="E [keV]", title="East drift")
scatter!(ax1, dx, E; markersize=6, label="Samples", color=:grey, alpha=0.5)
scatter!(ax1, dx_med, E_keV; color="#E69F00", label="Median")
lines!(ax1, dx_fit.(E_keV), E_keV; color=:black, linewidth=2)
axislegend(ax1; position=:lt)

ax2 = Axis(fig[1, 2]; xlabel="t [s]", ylabel="E [keV]", title="Arrival time for onset")
scatter!(ax2, t, E; markersize=6, label="Samples", color=:grey, alpha=0.5)

#scatter!(ax2, t_min, E_keV; color= "#56B4E9", label="minimum")
#lines!(ax2, tmin_fit.(E_keV), E_keV; color=:black, linewidth=2)

#scatter!(ax2, t_q01, E_keV; color="#56B4E9", label="1%")
#lines!(ax2, tq01_fit.(E_keV), E_keV; color=:black, linewidth=2)

scatter!(ax2, t_q05, E_keV; color= "#56B4E9", label="5%")
lines!(ax2, tq05_fit.(E_keV), E_keV; color=:black, linewidth=2)

axislegend(ax2; position=:lt)


linkyaxes!(ax1, ax2)
hideydecorations!(ax2, grid = false)

ylims!(ax1, 7, 40)
xlims!(ax1, 0, maximum(dx))
xlims!(ax2, 0, maximum(t))

#save("src/WPI/test_scripts/input_flux/figures/11262026_DriftAndTimeVsEnergy.png", fig)
