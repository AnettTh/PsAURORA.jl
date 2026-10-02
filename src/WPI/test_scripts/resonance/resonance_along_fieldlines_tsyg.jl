using AURORA
using CairoMakie


# NOTE: This does not work
## Make the data
L_vals = [3, 4, 5, 6, 7, 8]
ϕ = 0.0
colors = [:blue, :red, :green, :orange, :purple, :brown]

θ = deg2rad(0.0)
n = 1
E_eV = 30e3
μ = -cos(deg2rad(3.0))

ω_frac = 0.1

## Define figure
fig = Figure(size=(700, 400))
ax  = Axis(fig[1, 1];
    xlabel = "X [RE]",
    ylabel = "Z [RE]",
    title  = "Resonance velocity as a function of latitude",
    aspect = DataAspect()
)

## Make Earth
Φ = range(0, 2π, length=100)
poly!(ax, Point2f.(cos.(Φ), sin.(Φ)); color=:black, strokecolor=:black, strokewidth=1)


## Compute
field_line_data = map(L_vals) do L

    # Starting position at 4 MLT
    r0 = [RE*L * cos(ϕ), RE*L * sin(ϕ), 0.0]
    particle = ParticleState(E_eV, μ, r0, tsyganenko_field; relativistic=true)

    ## Trace field line from starting point
    xs_f, ys_f, zs_f, _ = trace_with_density(tsyganenko_field, DentonDensity(tsyganenko_field, L, ϕ), r0...; ds= RE*0.02)
    xs_b, ys_b, zs_b, _ = trace_with_density(tsyganenko_field, DentonDensity(tsyganenko_field, L, ϕ), r0...; ds=-RE*0.02)

    xs = [reverse(xs_b); xs_f[2:end]]
    ys = [reverse(ys_b); ys_f[2:end]]
    zs = [reverse(zs_b); zs_f[2:end]]

    # λ at each point along field line
    rs  = @. sqrt(xs^2 + ys^2 + zs^2)
    λs  = @. asin(zs / rs)
    Ls  = @. rs / (RE * cos(λs)^2)

    # λ_grid from field line trace
    λ_grid_L = λs

    ## Construct PlasmaState on traced field line
    plasma_L = PlasmaState(λ_grid_L, ϕ, DentonDensity(tsyganenko_field, L, ϕ), tsyganenko_field, Float64(L))

    ω = ω_frac * plasma_L.Ω_e[1]

    V_R_grid = map(enumerate(λ_grid_L)) do (i, λ)
        Ω_e   = plasma_L.Ω_e[i]
        ω_pe  = plasma_L.ω_pe[i]
        k     = dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)
        k_par = k * cos(θ)

        iszero(k_par) && return NaN
        V_R = ((ω + n * Ω_e / particle.γ) / k_par) / c₀
        (V_R > 1.0 || V_R < 0.0) && return NaN
        return V_R
    end

    # Project onto meridional plane for plotting
    r_perp = @. -sqrt(xs^2 + ys^2) / RE   # negative = nightside
    z_plot = zs ./ RE

    return (; x=r_perp, z=z_plot, V_R_grid)
end

## Make the lines
all_vals = filter(!isnan, vcat([vec(d.V_R_grid) for d in field_line_data]...))
clims_log = (log10(minimum(abs.(all_vals))), log10(maximum(abs.(all_vals))))

# Plot
for d in field_line_data
    lines!(ax, d.x,  d.z; color=log10.(abs.(d.V_R_grid)), colormap=:turbo, colorrange=clims_log, linewidth=3)
    lines!(ax, d.x, -d.z; color=log10.(abs.(d.V_R_grid)), colormap=:turbo, colorrange=clims_log, linewidth=3)
end

## Add text
parameters_text = """
    n_e = Denton-model
    n = $n
    E = $(E_eV/1e3) keV
    α = $(round(rad2deg(acos(abs(μ))), digits=1))°
    ω = $(ω_frac) Ωₑ
    θ = 0°
    """

text!(ax, 0.02, 0.98;
    text       = parameters_text,
    space      = :relative,
    align      = (:left, :top),
    fontsize   = 11
)

##
tick_vals     = [0.0, 10^(-0.25), 10^(-0.5), 10^(-0.75), 1.0]
tick_positions = log10.(tick_vals)
tick_labels   = ["0.0", "0.75", "0.5", "0.25", "1.0"]

Colorbar(fig[1, 2];
    colormap   = :turbo,
    limits     = clims_log,
    label      = "Resonance velocity [v∥/c₀]",
    ticks      = (tick_positions, tick_labels),
    tellheight = false
)

##
xlims!(0, -8.5)
ylims!(0, 4)

##
#save("src/WPI/test_scripts/resonance/res_lat_$(round(ω_frac, digits=2))_tsyg.png", fig)
