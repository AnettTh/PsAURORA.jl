using AURORA
using CairoMakie
using LinearAlgebra

## Grid in X-Z-plane
x_grid = range(-8*RE, -1*RE, length=50)
z_grid = range(-5*RE,  5*RE, length=50)

## Seed points for field line tracing (at equator, different x values)
x_seeds = range(-2*RE, -8*RE, length=7)

## Trace field lines by following B direction
function trace_fieldline(magnetic_field, x0, z0; ds=RE*0.05, n_steps=500)
    xs = [x0]
    zs = [z0]
    Bs = Float64[]

    x, z = x0, z0
    for _ in 1:n_steps
        # Evaluate and store
        B = magnetic_field(x, 0.0, z)
        Bmag = norm(B)
        push!(Bs, Bmag)
        Bx, _, Bz = B ./ Bmag

        # Update position
        x += ds * Bx
        z += ds * Bz
        push!(xs, x)
        push!(zs, z)

        # Stop inside or too far away from the Earth
        norm([x, z]) < RE && break
        norm([x, z]) > 10*RE && break
    end

    # Make sure lengths match
    push!(Bs, Bs[end])
    return xs ./ RE, zs ./ RE, Bs
end

## Collect field line data
fieldline_data = map(x_seeds) do x0
    # Trace both directions (forward and backwards) from equator
    xs_f, zs_f, Bs_f = trace_fieldline(tsyganenko_field, x0, 0.0; ds= RE*0.05)
    xs_b, zs_b, Bs_b = trace_fieldline(tsyganenko_field, x0, 0.0; ds=-RE*0.05)

    # Combine both directions, making sure all values are stored once
    xs = [reverse(xs_b); xs_f[2:end]]
    zs = [reverse(zs_b); zs_f[2:end]]
    Bs = [reverse(Bs_b); Bs_f[2:end]]
    return (; xs, zs, Bs)
end

## Define the color-range
all_Bs = vcat([log10.(d.Bs) for d in fieldline_data]...)
clims  = extrema(filter(isfinite, all_Bs))

## Plot
fig = Figure(size=(700, 400))
ax  = Axis(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Magnetic field using the Tsyganenko model",
    aspect = DataAspect()
)

# Make the Earth
ϕ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(ϕ), sin.(ϕ)); color=:black, strokecolor=:black, strokewidth=1)

# Make the field-lines
for d in fieldline_data
    lines!(ax, d.xs, d.zs;
        color      = log10.(d.Bs),
        colormap   = :turbo,
        colorrange = clims,
        linewidth  = 2
    )
end

Colorbar(fig[1, 2];
    colormap   = :turbo,
    limits     = clims,
    label      = "log₁₀(|B|) [T]",
    tellheight = false
)


## Define magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Seed points at equator for different L-shells
L_vals = range(2, 8, length=7)

## Trace field lines in 3D, project onto meridional plane
function trace_fieldline(
    magnetic_field::AbstractMagneticField,
    x0,
    y0,
    z0,
    ϕ_eq;
    ds=RE*0.05,
    n_steps=500
    )

    xs, ys, zs = [x0], [y0], [z0]
    Bs = Float64[]

    x, y, z = x0, y0, z0
    for _ in 1:n_steps

        # Evaluate magnetic field
        B = magnetic_field(Cartesian(x, y, z))
        Bmag = norm(B)

        # Store values
        push!(Bs, Bmag)
        iszero(Bmag) && break
        Bx, By, Bz = B ./ Bmag

        # Update position along field-line
        x += ds * Bx
        y += ds * By
        z += ds * Bz
        push!(xs, x); push!(ys, y); push!(zs, z)

        # Stop if crossing the Earth surface or moving outside valid range
        norm([x, y, z]) < RE       && break
        norm([x, y, z]) > 10 * RE  && break
    end
    push!(Bs, Bs[end])

    # Check if night-side or day-side
    r_sign = (180 < rad2deg(ϕ_eq) ≤ 360) ? -1 : 1

    # Project onto meridional plane
    r_perp = @. r_sign * sqrt(xs^2 + ys^2) / RE
    z_plot = zs ./ RE
    return r_perp, z_plot, Bs
end

## Collect field line data
fieldline_data = map(L_vals) do L
    x0 = L * RE * cos(ϕ_eq)
    y0 = L * RE * sin(ϕ_eq)
    z0 = 0.0

    rp_f, zp_f, Bs_f = trace_fieldline(magnetic_field, x0, y0, z0, ϕ_eq; ds= RE*0.05)
    rp_b, zp_b, Bs_b = trace_fieldline(magnetic_field, x0, y0, z0, ϕ_eq; ds=-RE*0.05)

    r_perp = [reverse(rp_b); rp_f[2:end]]
    z_plot = [reverse(zp_b); zp_f[2:end]]
    Bs = [reverse(Bs_b); Bs_f[2:end]]
    return (; r_perp, z_plot, Bs)
end

## Define the colorrange
all_Bs = vcat([log10.(d.Bs) for d in fieldline_data]...)
clims  = extrema(filter(isfinite, all_Bs))

## Make the plot
fig = Figure(size=(700, 500))
ax  = Axis(fig[1, 1];
    xlabel = "Rₚₑᵣₚ [RE]",
    ylabel = "Z [RE]",
    title  = "Tsyganenko field lines at $(longitude_to_MLT(ϕ_eq; degrees=false)) MLT",
    aspect = DataAspect()
)

## Add the Earth
θ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)

## Add the field-lines
for d in fieldline_data
    lines!(ax, d.r_perp, d.z_plot;
        color      = log10.(d.Bs),
        colormap   = :turbo,
        colorrange = clims,
        linewidth  = 2
    )
end

Colorbar(fig[1, 2];
    colormap   = :turbo,
    limits     = clims,
    label      = "log₁₀(|B|) [T]",
    tellheight = false
)

xlims!(ax, 0, -9)
