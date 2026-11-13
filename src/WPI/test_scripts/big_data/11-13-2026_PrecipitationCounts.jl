using AURORA
using JLD2
using CairoMakie

## Get the data, TOF with the shape [E, α, ω, Ψ₀]
data = load("src/WPI/data/11-13-2026/SmallTestSet.jld2")
#data = load("src/WPI/data/11-13-2026/energy-time-position.jld2")

E_grid = data["E_grid"] ./ 1e3
α_grid = rad2deg.(data["α_grid"])
ω_fracs = data["ω_frac"]
Ψ0_grid = data["Ψ0_grid"]

TOF = data["TOF"]

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
for k in Base.axes(logged_time, 3), iα in 1:nα, iE in 1:nE
    t = logged_time[iE, iα, k]
    check_valid(t) || continue
    it = min(searchsortedlast(t_edges, t), n_bins)    # the last edge goes in the last bin
    N[iE, iα, it] += 1
end

## Separate α and E into #(E, t) and #(α, t)
N_Et = dropdims(sum(N, dims=2), dims=2)
N_αt = dropdims(sum(N, dims=1), dims=1)

# Total count per (E, α) cell, independent of arrival time
N_Eα = dropdims(sum(N, dims=3), dims=3)

# Make common colorrange
all_counts = vcat(vec(N_Et), vec(N_αt))
colorrange = (0, max(1, maximum(all_counts)))

## Make figure
fig1 = Figure(size=(1200, 450))

Label(fig1[0, 1:3], "For all Ψ₀ and ω, how many particles precipitate?"; fontsize=18, font=:bold, tellwidth=false)

ax1 = Axis(fig1[1, 1]; xlabel="time [s]", ylabel="E [keV]", title="All pitch angles")
hm1 = heatmap!(
    ax1,
    t_centers,
    E_grid,
    N_Et';
    colormap=:turbo,
    colorrange=colorrange,
)

ax2 = Axis(fig1[1, 2]; xlabel="time [s]", ylabel="α [°]", title="All energies")
hm2 = heatmap!(
    ax2,
    t_centers,
    α_grid,
    N_αt';
    colormap=:turbo,
    colorrange=colorrange,
)

Colorbar(fig1[1, 3], hm1; label = "# precipitations")
display(fig1)

##
fig2 = Figure(size=(700, 500))
ax_total = Axis(fig2[1, 1];
    xlabel = "α [°]",
    ylabel = "E [keV]",
    title  = "Total number of precipitating particles",
)
hm_total = heatmap!(
    ax_total,
    α_grid,
    E_grid,
    N_Eα';
    colormap=:turbo,
    colorrange=(0, max(1, maximum(N_Eα))),
)
Colorbar(fig2[1, 2], hm_total; label = "# precipitations out of $(nω * nΨ) possible")
display(fig2)


## Make time-evolution, this was meant to be some sort of movie?
for it in 1:n_bins
    counts = N[:, :, it]

    fig3 = Figure(size=(600, 450))

    ax3 = Axis(
        fig3[1, 1];
        xlabel="α [°]",
        ylabel="E [keV]",
        title="How many particles precipitate? TOF: $(round(t_edges[it], digits=2)) s"
    )
    hm3 = heatmap!(
        ax3,
        α_grid,
        E_grid,
        counts';
        colormap=:turbo,
        colorrange=colorrange,
    )

    Colorbar(fig3[1, 2], hm3; label = "# precipitations")
    display(fig3)
    #save("src/WPI/test_scripts/big_data/figures/time/count_at_t$(round(t_centers[it], digits=2)).png", fig)
end
