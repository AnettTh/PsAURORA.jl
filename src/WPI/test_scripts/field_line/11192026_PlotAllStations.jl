using AURORA
using WGLMakie
using Bonito
using Dates

## Make magnetic field
t = DateTime("2020-12-21T03:00:00")
stations = [:tromso, :longyearbyen, :sodankyla]
station_colors = ["#E69F00", "#CC79A7", "#F0E442"]

t = DateTime("2020-12-21T03:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 4.0

magnetic_field = TsyganenkoMagneticField(
    t;
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)
## Field line through each station (south footpoint -> north footpoint, SM frame)
station_lines = []
for station in stations
    try
        push!(station_lines, (station=station, line=StationFieldLine(magnetic_field, LTPReference(station), t)))
    catch e
        @warn "No closed field line for $station" exception=e     # e.g. open or very long lines
    end
end

## Background field lines, seeded at the ionosphere on a grid of SM latitude and longitude
r_iono   = RE + z_ionosphere
seed_lat = deg2rad.(55.0:5.0:75.0)
seed_lon = deg2rad.(0.0:30.0:330.0)
seeds    = [(λ, ϕ) for λ in seed_lat, ϕ in seed_lon]
traces   = Vector{Any}(nothing, length(seeds))

ds_plot = 0.05 * RE                                  # coarse step is enough for plotting

Threads.@threads :greedy for k in eachindex(seeds)
    λ, ϕ = seeds[k]
    P0 = Cartesian(r_iono * cos(λ) * cos(ϕ), r_iono * cos(λ) * sin(ϕ), r_iono * sin(λ))
    try
        xs, ys, zs, _ = trace_fieldline(magnetic_field, P0; ds=ds_plot, n_steps=2_000, r_limit=20RE)
        traces[k] = (xs=xs, ys=ys, zs=zs, R_apex=maximum(@. sqrt(xs^2 + ys^2 + zs^2)))
    catch
        traces[k] = nothing                          # open line or field model failure: skip
    end
end

## Make figure
fig = Figure(size=(900, 700))
ax  = Axis3(fig[1, 1];
    xlabel = "X (SM) [RE]",
    ylabel = "Y (SM) [RE]",
    zlabel = "Z (SM) [RE]",
    title  = "Tsyganenko field lines at $(t)",
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

# Thin background field lines, coloured by the largest radius they reach
valid   = filter(!isnothing, traces)
R_range = extrema(tr.R_apex / RE for tr in valid)
for tr in valid
    pts = Point3f.(tr.xs ./ RE, tr.ys ./ RE, tr.zs ./ RE)
    lines!(ax, pts;
        color      = fill(tr.R_apex / RE, length(pts)),
        colormap   = :blues,
        colorrange = R_range,
        linewidth  = 1,
    )
end
Colorbar(fig[1, 2]; colormap=:blues, limits=R_range, label="Largest radius [RE]")

## Field lines through the station reference points
for (st, color) in zip(station_lines, station_colors)
    line = st.line
    lines!(ax,
        line.xs ./ RE,
        line.ys ./ RE,
        line.zs ./ RE;
        color = color,
        linewidth = 3,
        label = "R_max = $(round(line.r_max / RE, digits=2)) RE ($(st.station))",
    )

    # Magnetic equator (minimum B)
    scatter!(ax,
        [line.position_eq.x / RE],
        [line.position_eq.y / RE],
        [line.position_eq.z / RE];
        color = color, markersize = 10, marker = :circle,
    )

    # Ionosphere above the station: the last point, since the arrays run south to north
    scatter!(ax,
        [line.xs[end] / RE],
        [line.ys[end] / RE],
        [line.zs[end] / RE];
        color = color, markersize = 10, marker = :rect,
    )
end

axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
