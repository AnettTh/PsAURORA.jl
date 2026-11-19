using AURORA
using Dates
using WGLMakie
using Bonito
using LinearAlgebra

## Make dipole field
dipole = DipoleMagneticField()

## Set values for Tsyganenko field
times = [DateTime("2020-12-21T00:00:00") + Hour(h) for h in (1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)]
pdyn = 2.0
dst = -50
byimf = 1.0
bzimf = 2.0
station = :sodankyla

## Define wanted footpoint
ref = LTPReference(station)

function try_line(field, ref, t)
    try
        return StationFieldLine(field, ref, t)
    catch e
        @warn "No closed field line at $t" exception=e
        return nothing
    end
end

## Run function for each time
lines_tsyg   = []
lines_dipole = []
for t in times
    tsyg = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)
    push!(lines_tsyg,   try_line(tsyg,   ref, t))
    push!(lines_dipole, try_line(dipole, ref, t))
end


## Make figure
fig = Figure(size=(900, 700))
ax  = Axis3(fig[1,1];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    zlabel = "Z [RE]",
    aspect = :data,
    title = "Magnetic field at $(Date(first(times))) for $(station)"
)

## Add Earth sphere
u = range(0, π,  length=30)
v = range(0, 2π, length=30)
surface!(ax,
    [sin(ui)*cos(vi) for ui in u, vi in v],
    [sin(ui)*sin(vi) for ui in u, vi in v],
    [cos(ui)         for ui in u, vi in v];
    colormap = :blues, alpha = 0.3, shading = NoShading
)


## Divide MLT into colors
colors = cgrad(:turbo, length(times), categorical=true)

## Plot by looping through every trace
for (i, t) in enumerate(times)
    line_tsyg   = lines_tsyg[i]
    line_dipole = lines_dipole[i]

    # MLT of the station, read from the north footpoint (the last point; arrays run south to north)
    ref_line = something(line_tsyg, line_dipole, Some(nothing))
    mlt = isnothing(ref_line) ? NaN :
        longitude_to_MLT(atan(ref_line.ys[end], ref_line.xs[end]); degrees=true)
    mlt_str = round(mlt, digits=1)

    # Tsyganenko line
    if isnothing(line_tsyg)
        # Nothing to draw, but keep the legend entry
        lines!(ax, [NaN], [NaN], [NaN]; color=colors[i], linewidth=2,
               label="MLT $(mlt_str), open fieldline")
    else
        lines!(ax,
            line_tsyg.xs ./ RE, line_tsyg.ys ./ RE, line_tsyg.zs ./ RE;
            color=colors[i],
            linewidth=2,
            label="MLT $(mlt_str), R_max=$(round(line_tsyg.r_max / RE, digits=1)) RE",
        )
        # Magnetic equator (minimum B)
        scatter!(ax,
            [line_tsyg.position_eq.x / RE],
            [line_tsyg.position_eq.y / RE],
            [line_tsyg.position_eq.z / RE];
            color=colors[i], markersize=10, marker=:circle,
        )
    end

    # Dipole line
    if !isnothing(line_dipole)
        lines!(ax,
            line_dipole.xs ./ RE, line_dipole.ys ./ RE, line_dipole.zs ./ RE;
            color=colors[i],
            linewidth=1,
            linestyle=:dash,
        )
        scatter!(ax,
            [line_dipole.position_eq.x / RE],
            [line_dipole.position_eq.y / RE],
            [line_dipole.position_eq.z / RE];
            color=colors[i], markersize=5, marker=:circle,
        )
        # Start point: the station, at the north end of the line
        scatter!(ax,
            [line_dipole.xs[end] / RE],
            [line_dipole.ys[end] / RE],
            [line_dipole.zs[end] / RE];
            color=colors[i], markersize=12, marker=:circle,
        )
    end
end

lines!(ax, [NaN], [NaN], [NaN], color=:black, label="Tsyganenko")
lines!(ax, [NaN], [NaN], [NaN], color=:black, linestyle=:dash, label="Dipole")

axislegend(ax)
##
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
