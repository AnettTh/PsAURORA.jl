using AURORA
using CairoMakie
using LinearAlgebra

L = 6.5
ϕ_eq = deg2rad(240.0)

magnetic_field = TsyganenkoMagneticField(ϕ_eq)
#magnetic_field = DipoleMagneticField()
##
B_interp, s_grid = get_magnetic_field(magnetic_field, L, ϕ_eq)

## Compute B magnitude along field line
B_vals = [norm(B_interp(s)) for s in s_grid]

fig = Figure(size=(700, 400))
ax  = Axis(fig[1, 1];
    xlabel = "Arc length [RE]",
    ylabel = "|B| [T]",
    title  = "Magnetic field along field line (L=$(6.5))"
)

lines!(ax, s_grid ./ RE, B_vals)

# Mark equator (minimum B)
idx_eq = argmin(B_vals)
scatter!(ax, [s_grid[idx_eq] / RE], [B_vals[idx_eq]];
    color=:red, markersize=10, label="Equator")

axislegend(ax)
fig


##
fig = Figure(size=(800, 400))
ax  = Axis(fig[1,1]; xlabel="Arc length [RE]", ylabel="|B| [T]", title="L=6.5")

B_itp_dip,  s_dip  = get_magnetic_field(dipole,     L, ϕ_eq)
B_itp_tsyg, s_tsyg = get_magnetic_field(tsyganenko, L, ϕ_eq)

lines!(ax, s_dip  ./ RE, [norm(B_itp_dip(s))  for s in s_dip];  label="Dipole")
lines!(ax, s_tsyg ./ RE, [norm(B_itp_tsyg(s)) for s in s_tsyg]; label="Tsyganenko")
axislegend(ax)
fig


##
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Arc length [RE]",
    ylabel="Field strength [T]",
    title="L=$(L), MLT=$(round(longitude_to_MLT(ϕ_eq; degrees=false), digits=1))"
)

lines!(ax, s_tsyg ./ RE, [B_itp_tsyg(s)[1] for s in s_tsyg], label="Bx")
lines!(ax, s_tsyg ./ RE, [B_itp_tsyg(s)[2] for s in s_tsyg], label="By")
lines!(ax, s_tsyg ./ RE, [B_itp_tsyg(s)[3] for s in s_tsyg], label="Bz")

axislegend(ax; position=:ct)

fig
