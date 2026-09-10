## ===== Setup ===== ##
using CairoMakie
using AURORA
using AURORA: RE, qₑ, mₑ
using LinearAlgebra

L = 6.5
E_eV = 1e4
μ = 0.99

## ===== Field line ===== ##
function dipole_fieldline(L; n_points=500)
    λ_max = acos(sqrt(1 / L))
    λ = range(-λ_max, λ_max, length=n_points)
    r = L * RE .* cos.(λ).^2
    x = r .* cos.(λ)
    z = r .* sin.(λ)
    return x, z
end

x_fl, z_fl = dipole_fieldline(L)

## ===== Sample points along field line ===== ##
n_arrows = 10
λ_max = acos(sqrt(1 / L))
λ_samples = range(-λ_max * 0.9, λ_max * 0.9, length=n_arrows)

r_samples = L * RE .* cos.(λ_samples).^2
x_samples = r_samples .* cos.(λ_samples)
z_samples = r_samples .* sin.(λ_samples)

arrow_scale = RE * 0.3

## ===== Compute basis vectors ===== ##
bx_false = Float64[]; bz_false = Float64[]
bx_true  = Float64[]; bz_true  = Float64[]
vx_false = Float64[]; vz_false = Float64[]
vx_true  = Float64[]; vz_true  = Float64[]

for (xi, zi) in zip(x_samples, z_samples)
    r0_i = [xi, 0.0, zi]   # ← renamed to r0_i
    B = dipole_field(xi, 0.0, zi)

    b_f, _, _, _ = magnetic_basis(B; towards_equator=false)
    b_t, _, _, _ = magnetic_basis(B; towards_equator=true)
    # ...
    v0_f = get_v0_from_Eμ(dipole_field, r0_i, E_eV, μ; towards_equator=false)
    v0_t = get_v0_from_Eμ(dipole_field, r0_i, E_eV, μ; towards_equator=true)

    v_hat_f = collect(v0_f) ./ norm(collect(v0_f))
    v_hat_t = collect(v0_t) ./ norm(collect(v0_t))

    push!(vx_false, v_hat_f[1] * arrow_scale)
    push!(vz_false, v_hat_f[3] * arrow_scale)
    push!(vx_true,  v_hat_t[1] * arrow_scale)
    push!(vz_true,  v_hat_t[3] * arrow_scale)
end

## ===== Figure ===== ##
fig = Figure(resolution=(1600, 800))
θ = range(0, 2π, length=100)

# --- Panel 1: b̂ --- #
ax1 = Axis(fig[1, 1];
    xlabel="X [RE]", ylabel="Z [RE]",
    title="b̂ direction (L = $L)",
    aspect=DataAspect()
)
lines!(ax1, x_fl ./ RE, z_fl ./ RE, color=:lightgray, linewidth=2)
poly!(ax1, Point2f.(cos.(θ), sin.(θ)), color=:lightblue, strokecolor=:black, strokewidth=1)
arrows2d!(ax1, x_samples ./ RE, z_samples ./ RE, bx_false ./ RE, bz_false ./ RE,
    color=:blue, label="towards_equator=false")
arrows2d!(ax1, x_samples ./ RE, z_samples ./ RE, bx_true ./ RE, bz_true ./ RE,
    color=:red, label="towards_equator=true")
scatter!(ax1, x_samples ./ RE, z_samples ./ RE, color=:black, markersize=6)
axislegend(ax1, position=:rt)

# --- Panel 2: v0 --- #
ax2 = Axis(fig[1, 2];
    xlabel="X [RE]", ylabel="Z [RE]",
    title="v0 direction (μ = $μ, L = $L)",
    aspect=DataAspect()
)
lines!(ax2, x_fl ./ RE, z_fl ./ RE, color=:lightgray, linewidth=2)
poly!(ax2, Point2f.(cos.(θ), sin.(θ)), color=:lightblue, strokecolor=:black, strokewidth=1)
arrows2d!(ax2, x_samples ./ RE, z_samples ./ RE, vx_false ./ RE, vz_false ./ RE,
    color=:blue, label="towards_equator=false")
arrows2d!(ax2, x_samples ./ RE, z_samples ./ RE, vx_true ./ RE, vz_true ./ RE,
    color=:red, label="towards_equator=true")
scatter!(ax2, x_samples ./ RE, z_samples ./ RE, color=:black, markersize=6)
axislegend(ax2, position=:rt)

save("magnetic_basis_comparison.png", fig, px_per_unit=2)
println("Saved: magnetic_basis_comparison.png")
