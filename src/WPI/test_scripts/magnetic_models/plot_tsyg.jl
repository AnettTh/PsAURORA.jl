using AURORA
using CairoMakie
using LinearAlgebra

## Grid in x-z plane (y=0)
x_grid = range(-8*RE, -1*RE, length=50)   # nightside
z_grid = range(-5*RE,  5*RE, length=50)

## Seed points for field line tracing (at equator, different x values)
x_seeds = range(-2*RE, -8*RE, length=7)

## Trace field lines by following B direction
function trace_fieldline(B_func, x0, z0; ds=RE*0.05, n_steps=500)
    xs = [x0]
    zs = [z0]
    Bs = Float64[]

    x, z = x0, z0
    for _ in 1:n_steps
        B = B_func(x, 0.0, z)
        Bmag = norm(B)
        push!(Bs, Bmag)
        Bx, _, Bz = B ./ Bmag
        x += ds * Bx
        z += ds * Bz
        push!(xs, x)
        push!(zs, z)
        norm([x, z]) < RE && break   # stop inside Earth
        norm([x, z]) > 10*RE && break  # stop far away
    end
    push!(Bs, Bs[end])  # match length
    return xs ./ RE, zs ./ RE, Bs
end

## Collect field line data
fieldline_data = map(x_seeds) do x0
    # trace both directions from equator
    xs_f, zs_f, Bs_f = trace_fieldline(tsyganenko_field, x0, 0.0; ds= RE*0.05)
    xs_b, zs_b, Bs_b = trace_fieldline(tsyganenko_field, x0, 0.0; ds=-RE*0.05)

    xs = [reverse(xs_b); xs_f[2:end]]
    zs = [reverse(zs_b); zs_f[2:end]]
    Bs = [reverse(Bs_b); Bs_f[2:end]]
    return (; xs, zs, Bs)
end

## Colorrange
all_Bs = vcat([log10.(d.Bs) for d in fieldline_data]...)
clims  = extrema(filter(isfinite, all_Bs))

## Plot
fig = Figure(size=(700, 400))
ax  = Axis(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Magnetic field -- Tsyganenko model",
    aspect = DataAspect()
)

# Earth
ϕ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(ϕ), sin.(ϕ)); color=:black, strokecolor=:black, strokewidth=1)

# Field lines
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


# TODO: Verify this!
## Magnetic field model
ϕ_MLT  = deg2rad(120.0)   # 4 MLT
B_tsyg = TsyganenkoMagneticField(ϕ_MLT)

## Seed points at equator for different L-shells
L_vals = range(2, 8, length=7)

## Trace field lines in 3D, project onto meridional plane
function trace_fieldline(B_model::AbstractMagneticField, x0, y0, z0;
                          ds=RE*0.05, n_steps=500)
    xs, ys, zs = [x0], [y0], [z0]
    Bs = Float64[]

    x, y, z = x0, y0, z0
    for _ in 1:n_steps
        B    = B_model(Cartesian(x, y, z))
        Bmag = norm(B)
        push!(Bs, Bmag)
        iszero(Bmag) && break
        Bx, By, Bz = B ./ Bmag
        x += ds * Bx
        y += ds * By
        z += ds * Bz
        push!(xs, x); push!(ys, y); push!(zs, z)
        norm([x, y, z]) < RE       && break
        norm([x, y, z]) > 10 * RE  && break
    end
    push!(Bs, Bs[end])

    # Project onto meridional plane
    r_perp = @. -sqrt(xs^2 + ys^2) / RE   # negative = nightside
    z_plot = zs ./ RE
    return r_perp, z_plot, Bs
end

## Collect field line data
fieldline_data = map(L_vals) do L
    x0 = L * RE * cos(ϕ_MLT)
    y0 = L * RE * sin(ϕ_MLT)
    z0 = 0.0

    rp_f, zp_f, Bs_f = trace_fieldline(B_tsyg, x0, y0, z0; ds= RE*0.05)
    rp_b, zp_b, Bs_b = trace_fieldline(B_tsyg, x0, y0, z0; ds=-RE*0.05)

    r_perp = [reverse(rp_b); rp_f[2:end]]
    z_plot = [reverse(zp_b); zp_f[2:end]]
    Bs     = [reverse(Bs_b); Bs_f[2:end]]
    return (; r_perp, z_plot, Bs)
end

## Colorrange
all_Bs = vcat([log10.(d.Bs) for d in fieldline_data]...)
clims  = extrema(filter(isfinite, all_Bs))

## Plot
fig = Figure(size=(700, 500))
ax  = Axis(fig[1, 1];
    xlabel = "R_⊥ [RE]",
    ylabel = "Z [RE]",
    title  = "Tsyganenko field lines — 4 MLT",
    aspect = DataAspect()
)

# Earth
θ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(θ), sin.(θ)); color=:black, strokecolor=:black, strokewidth=1)

# Field lines
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
