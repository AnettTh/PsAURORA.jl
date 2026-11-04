using AURORA
using JLD2
using CairoMakie

## Get the data, TOF with the shape [E, α, ω, Ψ₀]
data = load("src/WPI/data/11-04-2026/SmallTestSet.jld2")

E_grid = data["E_grid"] ./ 1e3
α_grid = rad2deg.(data["α_grid"])
ω_fracs = data["ω_frac"]
Ψ0_grid = data["Ψ0_grid"]

## Get L=6.5 from the TOF matrix (E=10, α=2, ω=10, L=2, Ψ0=8)
TOF = data["TOF"][:, :, :, :]

## Merge ω and Ψ₀ into a "samples axis"
nE, nα, nω, nΨ = size(TOF)
logged_time = reshape(TOF, nE, nα, nω * nΨ)    # size [E, α, #]

## Check number of usable samples
function check_valid(t)
    return !isnan(t) && t > 0.0
end

precipitation_time = [t for t in logged_time if check_valid(t)]
println("$(length(precipitation_time)) precipitations out of $(length(logged_time)) possible")

## Make time-bins
n_bins = 5      # Figure out how to get this accurate
t_edges = range(minimum(precipitation_time), maximum(precipitation_time), length=n_bins + 1)
t_centers = collect(0.5 .* (t_edges[1:end-1] .+ t_edges[2:end]))

## Build histogram of arrivals by counting #(E, alpha, t)
N = zeros(nE, nα, n_bins)
for k in axes(logged_time, 3), iα in 1:nα, iE in 1:nE
    t = logged_time[iE, iα, k]
    check_valid(t) || continue
    it = min(searchsortedlast(t_edges, t), n_bins)    # the last edge goes in the last bin
    N[iE, iα, it] += 1
end

## Separate α and E into #(E, t) and #(α, t)
N_Et = dropdims(sum(N, dims=2), dims=2)
N_αt = dropdims(sum(N, dims=1), dims=1)

# Make common colorrange
all_counts = vcat(vec(N_Et), vec(N_αt))
colorrange = (minimum(all_counts), maximum(all_counts))

## Make figure
fig = Figure(size=(1200, 450))

Label(fig[0, 1:3], "For all Ψ₀ and ω, how many particles precipitate?"; fontsize=18, font=:bold, tellwidth=false)

ax1 = Axis(fig[1, 1]; xlabel="time [s]", ylabel="E [keV]", title="All pitch angles")
hm1 = heatmap!(
    ax1,
    t_centers,
    E_grid,
    N_Et';
    colormap=:turbo,
    colorrange=colorrange,
)

ax2 = Axis(fig[1, 2]; xlabel="time [s]", ylabel="α [°]", title="All energies")
hm2 = heatmap!(
    ax2,
    t_centers,
    α_grid,
    N_αt';
    colormap=:turbo,
    colorrange=colorrange,
)

Colorbar(fig[1, 3], hm1; label = "# precipitations")


## Make time-evolution
for it in 1:n_bins
    counts = N[:, :, it]
    colorrange = (1, maximum(N))

    fig = Figure(size=(600, 450))

    ax = Axis(fig[1, 1]; xlabel="α [°]", ylabel="E [keV]", title="How many particles precipitate?")
    hm = heatmap!(
        ax,
        α_grid,
        E_grid,
        counts';
        colormap=:turbo,
        colorrange=colorrange,
    )

    Colorbar(fig[1, 2], hm; label = "# precipitations")
    save("src/WPI/test_scripts/big_data/figures/time/count_at_t$(round(t_centers[it], digits=2)).png", fig)
end
