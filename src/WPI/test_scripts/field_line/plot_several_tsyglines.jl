using AURORA
using Dates
using WGLMakie
using Bonito
using LinearAlgebra

## Make dipole field
dipole = DipoleMagneticField()

## Set values for Tsyganenko field
times = [DateTime("2020-10-01T00:00:00") + Hour(h) for h in (1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)]
pdyn = 2.0
dst = -50
byimf = 1.0
bzimf = 2.0
station = :tromso

## Define wanted footpoint
ref = LTPReference(station)
#r_start = RE+10

## Run function for each time
traces_tsyg   = []
traces_dipole = []
for t in times
    tsyg = TsyganenkoMagneticField(ref; time=t, dst=dst, bzimf=bzimf, byimf=byimf)
    push!(traces_tsyg,   find_L_shell(tsyg,   ref, t; r_start=r_start))
    push!(traces_dipole, find_L_shell(dipole, ref, t; r_start=r_start))
end


## Make figure
fig = Figure(size=(900, 700))
ax  = Axis3(fig[1,1];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    zlabel = "Z [RE]",
    aspect = :data,
    title = "Magnetic field at $(Date(time)) for $(station)"
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

## Add dipole axis
#lines!(
#    ax,
#    [0.0, 0.0],
#    [0.0, 0.0],
#    [0.0, 5.0];
#    color=:black,
#    linewidth=1,
#    label="Dipole axis",
#)

## Divide MLT into colors
colors = cgrad(:turbo, length(times), categorical=true)

## Plot by looping through every trace
for (i, t) in enumerate(times)
    # Get the results
    xs_tsyg,   ys_tsyg,   zs_tsyg,   L_tsyg   = traces_tsyg[i]
    xs_dipole, ys_dipole, zs_dipole, L_dipole = traces_dipole[i]

    rs_tsyg = sqrt.(xs_tsyg.^2 .+ ys_tsyg.^2 .+ zs_tsyg.^2) ./ RE
    rs_dip  = sqrt.(xs_dipole.^2 .+ ys_dipole.^2 .+ zs_dipole.^2) ./ RE

    # MLT of the start point, read from the trace (the trace is in SM)
    mlt = mod(12 + rad2deg(atan(ys_tsyg[1], xs_tsyg[1])) / 15, 24)

    # Field line trace
    lines!(
        ax,
        xs_tsyg./RE,
        ys_tsyg./RE,
        zs_tsyg./RE;
        color=colors[i],
        linewidth=2,
        label = L_tsyg > 14.9 ? "MLT $(round(mlt, digits=1)), Open fieldline" : "MLT $(round(mlt, digits=1)), L=$(round(L_tsyg, digits=1))"
    )
    lines!(
        ax,
        xs_dipole./RE,
        ys_dipole./RE,
        zs_dipole./RE;
        color=colors[i],
        linewidth=1,
        linestyle=:dash,
    )

    # Maximum radius point (L-shell)
    idx_max_tsyg = argmax(rs_tsyg)
    scatter!(ax,
        [xs_tsyg[idx_max_tsyg]/RE],
        [ys_tsyg[idx_max_tsyg]/RE],
        [zs_tsyg[idx_max_tsyg]/RE];
        color=colors[i],
        markersize=10,
        marker=:circle,
    )
    idx_max_dip = argmax(rs_dip)
    scatter!(ax,
        [xs_dipole[idx_max_dip]/RE],
        [ys_dipole[idx_max_dip]/RE],
        [zs_dipole[idx_max_dip]/RE];
        color=colors[i],
        markersize=5,
        marker=:circle,
    )

    # Start point (Tromsø)
    scatter!(
        ax,
        [xs_dipole[1]/RE],
        [ys_dipole[1]/RE],
        [zs_dipole[1]/RE];
        color=colors[i],
        markersize=12,
        marker=:circle,
    )

    ## Add the sun direction
    #ps = dipole_tilt(t)
    #lines!(ax, [0.0, L_axis * cos(ps)], [0.0, 0.0], [0.0, L_axis * sin(ps)];
    #    color=colors[i], linewidth=1, linestyle=:dot,
    #    label= i == 1 ? "Sun direction (dotted)" : nothing,
    #)
end

lines!(ax, [NaN], [NaN], [NaN], color=:black, label="Tsyganenko")
lines!(ax, [NaN], [NaN], [NaN], color=:black, linestyle=:dash, label="Dipole")

axislegend(ax)
##
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
