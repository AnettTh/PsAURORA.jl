using CairoMakie
using AURORA
using LinearAlgebra

## ==================== Field line tracing ==================== ##
# Analytical dipole field line: r = L * RE * cos²(λ)
function dipole_fieldline(L; n_points=500)
    λ_max = acos(sqrt(1 / L))  # maximum latitude for this L-shell
    λ = range(-λ_max, λ_max, length=n_points)

    r = L * RE .* cos.(λ).^2
    x = r .* cos.(λ)
    z = r .* sin.(λ)
    return x, z
end


function trace_fieldline(x0, z0; ds=1e5, n_steps=1000)
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
    return xs, zs
end


L_shells = [2, 3, 4, 5, 6, 7, 8, 9, 10]

## ==================== Figure 1: Field lines ==================== ##
fig1 = Figure(resolution=(800, 800))
ax1 = Axis(fig1[1, 1];
    xlabel = "x (RE)",
    ylabel = "z (RE)",
    title  = "Dipole field lines",
    aspect = DataAspect()
)

for L in L_shells
    x, z = dipole_fieldline(L)
    lines!(ax1, x ./ RE, z ./ RE, color=:blue, linewidth=1)
end

# Trace numerical field lines from starting points on each L-shell
for L in L_shells
    # Start at equator for each L-shell
    x0 = L * RE
    z0 = 0.0

    # Trace northward (following B)
    xs, zs = trace_fieldline(x0, z0; ds=1e5, n_steps=2000)
    lines!(ax1, xs ./ RE, zs ./ RE, color=:red, linewidth=1, linestyle=:dash)

    # Trace southward (against B)
    xs, zs = trace_fieldline(x0, z0; ds=-1e5, n_steps=2000)
    lines!(ax1, xs ./ RE, zs ./ RE, color=:red, linewidth=1, linestyle=:dash)
end

# Add legend
axislegend(ax1,
    [LineElement(color=:blue), LineElement(color=:red, linestyle=:dash)],
    ["Analytical", "Numerical (dipole_field)"],
    position=:rt
)

# Draw Earth
θ = range(0, 2π, length=100)
poly!(ax1, cos.(θ), sin.(θ), color=:lightblue, strokecolor=:black, strokewidth=1)

save("src/WPI/diagnostic_figures/dipole_fieldlines.png", fig1)

## ==================== Figure 2: Field strength |B| ==================== ##
fig2 = Figure(resolution=(800, 800))
ax2 = Axis(fig2[1, 1];
    xlabel = "x (RE)",
    ylabel = "z (RE)",
    title  = "|B| (T)",
    aspect = DataAspect()
)

# Sample points along each field line
x_pts = Float64[]
z_pts = Float64[]
B_pts = Float64[]

for L in L_shells
    x, z = dipole_fieldline(L)
    for (xi, zi) in zip(x, z)
        r = sqrt(xi^2 + zi^2)
        r ≤ RE && continue
        r > 10 * RE && continue
        try
            B = dipole_field(xi, 0.0, zi)
            push!(x_pts, xi / RE)
            push!(z_pts, zi / RE)
            push!(B_pts, norm(B))
        catch
            continue
        end
    end
end

sc = scatter!(ax2, x_pts, z_pts;
    color      = B_pts,
    colormap   = :plasma,
    markersize = 3
)
Colorbar(fig2[1, 2], sc, label="|B| (T)")

# Draw Earth
poly!(ax2, cos.(θ), sin.(θ), color=:lightblue, strokecolor=:black, strokewidth=1)

save("src/WPI/diagnostic_figures/dipole_field_strength.png", fig2)

## ==================== Figure 3: Field vectors ==================== ##
fig3 = Figure(resolution=(800, 800))
ax3 = Axis(fig3[1, 1];
    xlabel = "x (RE)",
    ylabel = "z (RE)",
    title  = "Dipole field vectors",
    aspect = DataAspect()
)

# Sample fewer points for arrows
x_arr = Float64[]
z_arr = Float64[]
Bx_arr = Float64[]
Bz_arr = Float64[]

for L in L_shells
    x, z = dipole_fieldline(L; n_points=30)  # fewer points for arrows
    for (xi, zi) in zip(x, z)
        r = sqrt(xi^2 + zi^2)
        r ≤ RE && continue
        r > 10 * RE && continue
        try
            B = dipole_field(xi, 0.0, zi)
            b = B / norm(B)   # unit vector
            push!(x_arr, xi / RE)
            push!(z_arr, zi / RE)
            push!(Bx_arr, b[1])
            push!(Bz_arr, b[3])
        catch
            continue
        end
    end
end

arrows!(ax3, x_arr, z_arr, Bx_arr, Bz_arr;
    arrowsize  = 8,
    lengthscale = 0.15,
    arrowcolor  = :red,
    linecolor   = :red
)

# Draw field lines for reference
for L in L_shells
    x, z = dipole_fieldline(L)
    lines!(ax3, x ./ RE, z ./ RE, color=:lightgray, linewidth=1)
end

# Draw Earth
poly!(ax3, cos.(θ), sin.(θ), color=:lightblue, strokecolor=:black, strokewidth=1)

save("src/WPI/diagnostic_figures/dipole_field_vectors.png", fig3)

println("Saved: dipole_fieldlines.png, dipole_field_strength.png, dipole_field_vectors.png")
