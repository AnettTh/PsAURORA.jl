using AURORA
using WGLMakie
using Bonito
using LinearAlgebra
using Dates

## Make magnetic field
t = DateTime("2020-01-01T00:01:30")
station = :tromso
ref = LTPReference(station)
tsyg = TsyganenkoMagneticField(ref; time=t)

## Trace the field line from the reference point
xs_tsyg, ys_tsyg, zs_tsyg, L_tsyg = find_L_shell(tsyg, ref, t)
rs_tsyg = sqrt.(xs_tsyg.^2 .+ ys_tsyg.^2 .+ zs_tsyg.^2) ./ RE

## Background field lines, seeded in the equatorial plane
seed_L   = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
seed_ϕ   = deg2rad.(0.0:30.0:330.0)
seeds    = [(L, ϕ) for L in seed_L, ϕ in seed_ϕ]      # matrix of (L, ϕ) pairs
traces   = Vector{Any}(undef, length(seeds))

# Coarse step is enough for plotting; each trace is traced both ways from the seed
ds_plot = 0.05 * RE

Threads.@threads :greedy for k in eachindex(seeds)
    L, ϕ = seeds[k]
    seed = Cartesian(L * RE * cos(ϕ), L * RE * sin(ϕ), 0.0)
    try
        _, trace = find_R_max(tsyg, seed; ds=ds_plot, store_trace=true)
        traces[k] = (L=L, trace=trace)
    catch
        traces[k] = nothing        # skip seeds where the field model fails
    end
end

## Make figure
fig = Figure(size=(900, 700))
ax  = Axis3(fig[1,1];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    zlabel = "Z [RE]",
    title  = "Tsyganenko field lines",
    aspect = :data,
)

# Earth sphere
u = range(0, π,  length=30)
v = range(0, 2π, length=30)
surface!(ax,
    [sin(ui)*cos(vi) for ui in u, vi in v],
    [sin(ui)*sin(vi) for ui in u, vi in v],
    [cos(ui)         for ui in u, vi in v];
    colormap = :blues, alpha = 0.3, shading = NoShading,
)

# Thin background field lines, coloured by equatorial seed L
L_range = (minimum(seed_L), maximum(seed_L))
for t in traces
    isnothing(t) && continue
    pts = Point3f.(t.trace.xs ./ RE, t.trace.ys ./ RE, t.trace.zs ./ RE)
    lines!(ax, pts;
        color      = fill(t.L, length(pts)),
        colormap   = :viridis,
        colorrange = L_range,
        linewidth  = 1,
    )
end
Colorbar(fig[1, 2]; colormap=:viridis, limits=L_range, label="L(Z=0)")

# Field line through the Tromsø reference point
lines!(ax,
    xs_tsyg ./ RE,
    ys_tsyg ./ RE,
    zs_tsyg ./ RE;
    color = :black,
    linewidth = 3,
    label = "Fieldline from $(station)",
)

# Maximum radius point (L-shell)
idx_max_tsyg = argmax(rs_tsyg)
scatter!(ax,
    [xs_tsyg[idx_max_tsyg] / RE],
    [ys_tsyg[idx_max_tsyg] / RE],
    [zs_tsyg[idx_max_tsyg] / RE];
    color = :black, markersize = 10, marker = :circle,
    label = "Lₘₐₓ=$(round(L_tsyg, digits=2))",
)

# Start point (ionosphere above Tromsø)
scatter!(ax,
    [xs_tsyg[1] / RE],
    [ys_tsyg[1] / RE],
    [zs_tsyg[1] / RE];
    color = :black, markersize = 10, marker = :circle,
)

axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
