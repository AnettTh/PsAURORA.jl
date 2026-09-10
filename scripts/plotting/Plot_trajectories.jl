using AURORA
using CairoMakie

##
z_end = 600e3   # [m]
r = RE + z_end  # [m]

L = 6.0

λ = acos(sqrt(r / (L * RE)))

r0 = [
    r * cos(λ),
    0.0,
    r * sin(λ)
]

r_source = L * RE

E_eV = 1e4

# θ = 180° (field-aligned towards Earth) means μ = -1
μ = -0.99

v0 = get_v0_from_Eμ(dipole_field, r0, E_eV, μ; ϕ=0.0, towards_equator=true)

##
result = boris_mover_TOF(
    dipole_field,
    r0,
    v0,
    r_source;
    n_T=100_000,
    resolution=10,
    store_trajectory=true
)


##
pos = result.r


##
fig = Figure(size=(800,400))

ax1 = Axis(
    fig[1, 1],
    xlabel="X [km]",
    ylabel="Z [km]"
)

lines!(
    ax1,
    pos[2:32, 1] ./ 1e3,
    pos[2:32, 3] ./ 1e3,
    label="Start",
    color=:red
)

scatter!(
    ax1,
    [pos[2, 1] / 1e3],
    [pos[2, 3] / 1e3],
    label="r0",
    color=:green,
    markersize=15
)

ax2 = Axis(
    fig[1, 2],
    xlabel="X [km]",
    ylabel="Z [km]"
)

lines!(
    ax2,
    pos[end-30:end, 1] ./ 1e3,
    pos[end-30:end, 3] ./ 1e3,
    label="End",
    color=:blue
)

fig[0, 1:2] = Legend(fig, [ax1, ax2], merge=true, orientation=:horizontal)

fig
