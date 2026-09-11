using CairoMakie
using AURORA
using LinearAlgebra

## Define parameters
L = 6.5
h_top = 600e3

r_top = RE + h_top
λ = acos(sqrt(r_top / (L * RE)))
r0 = [-r_top * cos(λ), 0.0, r_top * sin(λ)]

E_eV = 1e4
μ = 0.5

## Compute basis and velocity
B = dipole_field(r0...)
b, e1, e2, _ = magnetic_basis(B)

v0 = get_v0_from_Eμ(dipole_field, r0, E_eV, μ; towards_equator=true)
v0_hat = collect(v0) ./ norm(collect(v0))

## Make a field-line
function trace_fieldline(x0, z0; ds=1e4, n_steps=2000)
    x, z = x0, z0
    xs = [x]
    zs = [z]
    for _ in 1:n_steps
        B = dipole_field(x, 0.0, z)
        b = B / norm(B)
        x += ds * b[1]
        z += ds * b[3]
        r = sqrt(x^2 + z^2)
        r ≤ RE && break
        r > 10 * RE && break
        push!(xs, x)
        push!(zs, z)
    end
    return xs ./ RE, zs ./ RE
end


x_eq, z_eq = trace_fieldline(r0[1], r0[3]; ds=-1e4, n_steps=20000)
x_earth, z_earth = trace_fieldline(r0[1], r0[3]; ds=1e4, n_steps=20000)


## Make plot
s = RE * 0.4
p = RE * 1.2

fig = Figure(size=(800, 400))
ax = Axis(
    fig[1,1],
    xlabel="X [RE]",
    ylabel="Z [RE]",
    aspect=DataAspect(),
    xticks=-8:1:0,
    title="Magnetic basis and velocity")

ax.xreversed=true

## Make the Earth
θ = range(0, 2π, length=100)
lines!(ax, cos.(θ), sin.(θ), color=:black)
lines!(ax, (r_top/RE) * cos.(θ), (r_top/RE) * sin.(θ), color=:black, alpha=0.5)


## Make field-line
lines!(ax, x_eq, z_eq, color=:blue, linewidth=1)
lines!(ax, x_earth, z_earth, color=:blue, linewidth=1, label="Magnetic field")


## Magnetic basis
# b̂ — red
arrows2d!(ax, [r0[1]/RE], [r0[3]/RE], [b[1]*p/RE], [b[3]*p/RE],
    color=:red, label="b̂")

# e1 — blue
arrows2d!(ax, [r0[1]/RE], [r0[3]/RE], [e2[1]*p/RE], [e2[3]*p/RE],
    color=:blue, label="e₂")

# v0 — green
arrows2d!(ax, [r0[1]/RE], [r0[3]/RE], [v0_hat[1]*p/RE], [v0_hat[3]*p/RE],
    color=:green, label="v₀")

## Initial position
scatter!(ax, [r0[1]/RE], [r0[3]/RE], color=:black, markersize=8, label="r₀")

## Adjust what is shown on plot
ylims!(ax, 0, 3)
xlims!(ax, 1, -8)

## Legend and show
Legend(fig[1, 1],
    [
        LineElement(color=:red, linewidth=2),
        LineElement(color=:blue, linewidth=2),
        LineElement(color=:green, linewidth=2),
        MarkerElement(color=:black, marker=:circle, markersize=8)
    ],
    ["b̂", "e₂", "v₀", "r₀"],
    tellwidth=false,
    tellheight=false
)

## Save figure
save("magnetic_basis.png", fig)
