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

##
xlims!(ax, 0, -7)
