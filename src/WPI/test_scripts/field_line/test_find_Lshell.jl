using AURORA
using WGLMakie
using Bonito


## Make magnetic field
ϕ_eq = deg2rad(240.0)
tsyg = TsyganenkoMagneticField(ϕ_eq)
dipole = DipoleMagneticField()

## Define point of reference
#lat_ref = deg2rad(69.66)
#lon_ref = deg2rad(0.0)

##
ref = LTPReference(:tromso)
MAG = Geograpic2Magnetic(ref)

## Run function
xs_tsyg, ys_tsyg, zs_tsyg, L_tsyg = find_L_shell(
    tsyg,
    deg2rad(MAG.mlat_deg),
    deg2rad(MAG.mlon_deg),
)
xs_dipole, ys_dipole, zs_dipole, L_dipole = find_L_shell(
    dipole,
    deg2rad(MAG.mlat_deg),
    deg2rad(MAG.mlon_deg),
)

rs_tsyg = sqrt.(xs_tsyg.^2 .+ ys_tsyg.^2 .+ zs_tsyg.^2) ./ RE
rs_dip  = sqrt.(xs_dipole.^2  .+ ys_dipole.^2  .+ zs_dipole.^2)  ./ RE


##
fig = Figure(size=(800, 600))
ax  = Axis3(fig[1,1];
    xlabel = "X [RE]",
    ylabel = "Y [RE]",
    zlabel = "Z [RE]",
    title  = "Field line trace from Tromsø",
    aspect = :data
)

# Earth sphere
u = range(0, π,  length=30)
v = range(0, 2π, length=30)
surface!(ax,
    [sin(ui)*cos(vi) for ui in u, vi in v],
    [sin(ui)*sin(vi) for ui in u, vi in v],
    [cos(ui)         for ui in u, vi in v];
    colormap = :blues, alpha = 0.3, shading = NoShading
)

# Field line trace
lines!(
    ax,
    xs_tsyg./RE,
    ys_tsyg./RE,
    zs_tsyg./RE;
    color=:red,
    linewidth=2,
    label="Tsyganenko"
)
lines!(
    ax,
    xs_dipole./RE,
    ys_dipole./RE,
    zs_dipole./RE;
    color=:green,
    linewidth=2,
    label="Dipole"
)

# Maximum radius point (L-shell)
idx_max_tsyg = argmax(rs_tsyg)
scatter!(ax,
    [xs_tsyg[idx_max_tsyg]/RE],
    [ys_tsyg[idx_max_tsyg]/RE],
    [zs_tsyg[idx_max_tsyg]/RE];
    color=:red, markersize=10, marker=:circle,
    label="L=$(round(L_tsyg, digits=2)) (Tsyg)"
)
idx_max_dip = argmax(rs_dip)
scatter!(ax,
    [xs_dipole[idx_max_dip]/RE],
    [ys_dipole[idx_max_dip]/RE],
    [zs_dipole[idx_max_dip]/RE];
    color=:green, markersize=10, marker=:circle,
    label="L=$(round(L_dipole, digits=2)) (Dipole)"
)

# Start point (Tromsø)
scatter!(ax, [xs_dipole[1]/RE], [ys_dipole[2]/RE], [zs_dipole[2]/RE];
    color=:blue, markersize=15, marker=:circle, label="Ionosphere above Tromsø"
)


axislegend(ax)

isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
