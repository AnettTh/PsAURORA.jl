using AURORA
using CairoMakie
using LinearAlgebra: norm, dot

##
h_top = 600e3
r = RE + h_top

L = 6.0

λ = acos(sqrt(r / (L * RE)))

r0 = [
    r * cos(λ),
    0.0,
    r * sin(λ)
]

r_source = L * RE

E_eV = 1e4

μ = 0.99

v0 = get_v0_from_Eμ(dipole_field, r0, E_eV, μ; ϕ=0.0)

##
result = boris_mover_TOF(
    dipole_field,
    r0,
    v0,
    r_source;
    n_T=100_000,
    resolution=10
)


##
pos = result.r


##
fig = Figure()
ax = Axis(fig[1, 1])

lines!(
    ax,
    pos[2:end, 1] ./ RE,
    pos[2:end, 3] ./ RE,
)

scatter!(
    ax,
    [pos[2, 1] / RE],
    [pos[2, 3] / RE],
)

fig
