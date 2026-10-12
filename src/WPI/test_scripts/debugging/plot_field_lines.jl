using AURORA
using WGLMakie
using Bonito
using LinearAlgebra

## Make particle
L = 6.5
ϕ_eq = deg2rad(240.0)
μ = -cos(deg2rad(2.0))
E = 30e3
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))

p_dipole = ParticleState(
    E,
    μ,
    r0,
    DipoleMagneticField(),
    ϕ_eq;
    precomputed_magnetic_field=false,
    relativistic=true,
    φ0=deg2rad(0.0)
)

p_tsyg = ParticleState(
    E,
    μ,
    r0,
    TsyganenkoMagneticField(ϕ_eq),
    ϕ_eq;
    precomputed_magnetic_field=false,
    relativistic=true,
    φ0=deg2rad(0.0)
)

## Find trace
R_max_dipole, trace_dipole = find_R_max(p_dipole; store_trace=true)
R_max_tsyg, trace_tsyg = find_R_max(p_tsyg; store_trace=true)


## Compute resonance
λ_grid = range(-deg2rad(50), deg2rad(50), length=500)
plasma_dipole = PlasmaState(
    λ_grid,
    ϕ_eq,
    DentonDensity(DipoleMagneticField(), L, ϕ_eq),
    DipoleMagneticField(),
    Float64(L)
)
plasma_tsyg = PlasmaState(
    λ_grid,
    ϕ_eq,
    DentonDensity(TsyganenkoMagneticField(ϕ_eq), L, ϕ_eq),
    TsyganenkoMagneticField(ϕ_eq),
    Float64(L)
)

Ω_e_eq_dip  = minimum(plasma_dipole.Ω_e)
Ω_e_eq_tsyg = minimum(plasma_tsyg.Ω_e)

ω_dip  = Ω_e_eq_dip  * 0.3
ω_tsyg = Ω_e_eq_tsyg * 0.3

λ_res_dip,  V_r_dip  = resonance_latitude(ω_dip,  p_dipole, plasma_dipole; n=1)
λ_res_tsyg, V_r_tsyg = resonance_latitude(ω_tsyg, p_tsyg,   plasma_tsyg;   n=1)

## Helper to convert spherical to Cartesian RE
function λ_to_xyz(p, λ)
    r = Cartesian(Spherical(p.L, λ, p.ϕ_eq))
    return r.x/RE, r.y/RE, r.z/RE
end

## Boris mover trajectory
result_dip  = boris_mover(p_dipole; store_trajectory=true, n_T=100_000, λ_start=λ_res_dip)
result_tsyg = boris_mover(p_tsyg;   store_trajectory=true, n_T=100_000, λ_start=λ_res_tsyg)


## Make figure
fig = Figure()
ax = Axis3(
    fig[1,1],
    aspect=:equal,
    xlabel="X [RE]",
    ylabel="Y [RE]",
    zlabel="Z [RE]",
    title="Trace of magnetic field lines" *
        "L=$(L) at $(round(longitude_to_MLT(ϕ_eq; degrees=false), digits=1)) MLT"
)

# Starting positions
x0, y0, z0 = p_dipole.r0 ./ RE
scatter!(ax, [x0], [y0], [z0]; color=:red, markersize=20, marker=:star5, label="r0 (dipole)")
x0, y0, z0 = p_tsyg.r0 ./ RE
scatter!(ax, [x0], [y0], [z0]; color=:green, markersize=20, marker=:star5, label="r0 (tsyg)")

# Resonance latitude
if !isnan(λ_res_dip)
    x, y, z = λ_to_xyz(p_dipole, λ_res_dip)
    scatter!(ax, [x], [y], [z]; color=:red, markersize=20, marker=:diamond, label="λ_res (dipole)")
end
if !isnan(λ_res_tsyg)
    x, y, z = λ_to_xyz(p_tsyg, λ_res_tsyg)
    scatter!(ax, [x], [y], [z]; color=:green, markersize=20, marker=:diamond, label="λ_res (tsyg)")
end

# Ionospheric footpoint
for (p, color) in [(p_dipole, :red), (p_tsyg, :green)]
    for sign in [1, -1]
        x, y, z = λ_to_xyz(p, sign * p.λ_ionosphere)
        scatter!(ax, [x], [y], [z]; color=color, markersize=15, marker=:circle)
    end
end

# Boris mover trajectory
lines!(ax,
    result_dip.position[:, 1] ./ RE,
    result_dip.position[:, 2] ./ RE,
    result_dip.position[:, 3] ./ RE;
    color=:red, linewidth=1, linestyle=:dot, label="Boris (dipole)"
)
lines!(ax,
    result_tsyg.position[:, 1] ./ RE,
    result_tsyg.position[:, 2] ./ RE,
    result_tsyg.position[:, 3] ./ RE;
    color=:green, linewidth=1, linestyle=:dot, label="Boris (tsyg)"
)

# Magnetic equator
for (p, trace, color) in [(p_dipole, trace_dipole, :red), (p_tsyg, trace_tsyg, :green)]
    Bs   = [norm(p.magnetic_field(Cartesian(x, y, z)))
            for (x, y, z) in zip(trace.xs, trace.ys, trace.zs)]
    idx  = argmin(Bs)
    x_eq = trace.xs[idx] / RE
    y_eq = trace.ys[idx] / RE
    z_eq = trace.zs[idx] / RE
    scatter!(ax, [x_eq], [y_eq], [z_eq]; color=color, markersize=20, marker=:cross, label="B_eq ($(color))")
end

# Field-lines
lines!(
    ax,
    trace_dipole.xs ./ RE,
    trace_dipole.ys ./ RE,
    trace_dipole.zs ./ RE,
    label="Dipole",
    color=:red,
    alpha=0.5
)

lines!(
    ax,
    trace_tsyg.xs ./ RE,
    trace_tsyg.ys ./ RE,
    trace_tsyg.zs ./ RE,
    label="Tsyganenko",
    color=:green,
    alpha=0.5
)

# Earth sphere
u = range(0, π,  length=50)
v = range(0, 2π, length=50)
xs_sphere = [sin(ui)*cos(vi) for ui in u, vi in v]
ys_sphere = [sin(ui)*sin(vi) for ui in u, vi in v]
zs_sphere = [cos(ui)         for ui in u, vi in v]

surface!(ax, xs_sphere, ys_sphere, zs_sphere;
    colormap = :blues,
    alpha    = 0.5,
    shading  = NoShading
)

# Find R_max position from trace
rs_tsyg = sqrt.(trace_tsyg.xs.^2 .+ trace_tsyg.ys.^2 .+ trace_tsyg.zs.^2)
idx_max  = argmax(rs_tsyg)

x_max_tsyg = trace_tsyg.xs[idx_max] / RE
y_max_tsyg = trace_tsyg.ys[idx_max] / RE
z_max_tsyg = trace_tsyg.zs[idx_max] / RE

lines!(
    ax,
    [0.0, x_max_tsyg],
    [0.0, y_max_tsyg],
    [0.0, z_max_tsyg];
    color=:green,
    linewidth=2
)
scatter!(
    ax,
    [x_max_tsyg],
    [y_max_tsyg],
    [z_max_tsyg];
    color=:green,
    markersize=15
)

# Dipole R_max is always at equator
x_max_dip = R_max_dipole * cos(ϕ_eq) / RE
y_max_dip = R_max_dipole * sin(ϕ_eq) / RE
z_max_dip = 0.0

lines!(ax, [0.0, x_max_dip], [0.0, y_max_dip], [0.0, z_max_dip]; color=:red, linewidth=2)
scatter!(ax, [x_max_dip], [y_max_dip], [z_max_dip]; color=:red, markersize=15)


Legend(fig[1,2], ax)

##
lim = 8.0   # RE
xlims!(ax, -lim, 2)
ylims!(ax, -lim, 2)
zlims!(ax, -lim+3, lim-3)
## Refreshing the plot
#close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
