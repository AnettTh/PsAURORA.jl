using AURORA
using WGLMakie
using Bonito
using Dates

## Make magnetic field
t = DateTime("2020-12-21T03:00:00")
stations = [:tromso, :longyearbyen, :sodankyla]
station_colors = ["#E69F00", "#CC79A7", "#F0E442"]
tsyg = TsyganenkoMagneticField(LTPReference(stations[1]); time=t, r_limit=20RE)

## Trace the field line from each reference point
station_traces = []
for station in stations
    ref = LTPReference(station)
    res = find_L_shell(tsyg, ref, t)
    push!(station_traces, (station=station, xs=res.xs, ys=res.ys, zs=res.zs, L=res.L))
end

## Background field lines, seeded in the equatorial plane
seed_L   = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0]
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
for tr in traces
    isnothing(tr) && continue
    pts = Point3f.(tr.trace.xs ./ RE, tr.trace.ys ./ RE, tr.trace.zs ./ RE)
    lines!(ax, pts;
        color      = fill(tr.L, length(pts)),
        colormap   = :blues,
        colorrange = L_range,
        linewidth  = 1,
    )
end
Colorbar(fig[1, 2]; colormap=:blues, limits=L_range, label="L(Z=0)")

## Field lines through the station reference points
for (st, color) in zip(station_traces, station_colors)
    rs = sqrt.(st.xs.^2 .+ st.ys.^2 .+ st.zs.^2) ./ RE

    lines!(ax,
        st.xs ./ RE,
        st.ys ./ RE,
        st.zs ./ RE;
        color = color,
        linewidth = 3,
        label = "Lₘₐₓ=$(round(st.L, digits=2)) ($(st.station))",
    )

    # Maximum radius point (L-shell)
    idx_max = argmax(rs)
    scatter!(ax,
        [st.xs[idx_max] / RE],
        [st.ys[idx_max] / RE],
        [st.zs[idx_max] / RE];
        color = color, markersize = 10, marker = :circle,
    )

    # Start point (ionosphere above the station)
    scatter!(ax,
        [st.xs[1] / RE],
        [st.ys[1] / RE],
        [st.zs[1] / RE];
        color = color, markersize = 10, marker = :circle,
    )
end

axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
