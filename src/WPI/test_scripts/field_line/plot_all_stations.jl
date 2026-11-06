## Make magnetic field
t = DateTime("2020-01-01T00:01:30")
stations = [:tromso, :longyearbyen, :sodankyla]          # CHANGED: was a single station
station_colors = [:black, :red, :magenta]                # ADDED
tsyg = TsyganenkoMagneticField(LTPReference(stations[1]); time=t)    # CHANGED: see note 1

## Trace the field line from each reference point
station_traces = []                                      # ADDED
for station in stations                                  # ADDED
    ref = LTPReference(station)
    xs, ys, zs, L = find_L_shell(tsyg, ref, t)
    push!(station_traces, (station=station, xs=xs, ys=ys, zs=zs, L=L))
end

## Background field lines, seeded in the equatorial plane
# ... unchanged ...

## Make figure
# ... figure, axis, Earth sphere and background field lines unchanged ...
Colorbar(fig[1, 2]; colormap=:viridis, limits=L_range, label="L(Z=0)")

# Field lines through the station reference points
for (st, color) in zip(station_traces, station_colors)   # CHANGED: loop over stations
    rs = sqrt.(st.xs.^2 .+ st.ys.^2 .+ st.zs.^2) ./ RE

    lines!(ax,
        st.xs ./ RE,
        st.ys ./ RE,
        st.zs ./ RE;
        color = color,
        linewidth = 3,
        label = "Fieldline from $(st.station)",
    )

    # Maximum radius point (L-shell)
    idx_max = argmax(rs)
    scatter!(ax,
        [st.xs[idx_max] / RE],
        [st.ys[idx_max] / RE],
        [st.zs[idx_max] / RE];
        color = color, markersize = 10, marker = :circle,
        label = "Lₘₐₓ=$(round(st.L, digits=2)) ($(st.station))",
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
