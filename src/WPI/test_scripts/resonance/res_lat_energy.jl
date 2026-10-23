using AURORA
using CairoMakie

## Make magnetic field
ϕ = deg2rad(240.0)
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Parameters
L_vals = [4.0, 5.0, 6.0, 7.0]
μ      = -cos(deg2rad(3.0))
ω_frac = 0.35
n      = 1
θ      = deg2rad(0)
E_grid = range(1e3, 40e3, length=50)
colors = [:blue, :red, :green, :orange, :purple]

## Compute
fig = Figure(size=(600, 400))
ax  = Axis(fig[1, 1];
    xlabel = "E [keV]",
    ylabel = "λᵣₑₛ [°]",
    title  = "Resonance latitude vs energy \n" *
    if sign(μ) < 0 "propagating opposite to the wave" else "propagating in the same direction as the wave" end,
    #xscale = log10
)

##
for (L, color) in zip(L_vals, colors)
    @info "Starting L=$L"
    r0 = Cartesian(L*RE, 0.0, ϕ_eq)

    ## Define the plasma
    @info "Making plasma..."
    λ_grid   = range(0.0, deg2rad(60), length=500)
    plasma_L = PlasmaState(λ_grid, ϕ_eq, DentonDensity(magnetic_field, L, ϕ_eq), magnetic_field, Float64(L))

    ## Define the wave
    ω = ω_frac * minimum(plasma_L.Ω_e)

    @info "Computing resonance..."
    λ_res_grid = @showprogress map(E_grid) do E
        p = ParticleState(E, μ, r0, magnetic_field, ϕ_eq; relativistic=true)
        λ_res, Vr_res = resonance_latitude(ω, p, plasma_L; n=n, θ=θ)
        isnan(λ_res) ? NaN : rad2deg(λ_res)
    end

    lines!(ax, collect(E_grid)./1e3, λ_res_grid;
        color=color, label="L = $L")
end

## Add parameter text
text!(ax, 0.02, 0.98;
    text      = "ω = $(ω_frac)Ωe \n" *
                "α = $(round(rad2deg(acos(abs(μ))), digits=1))° \n" *
                "n = $n \n" *
                "θ = $(round(rad2deg(θ), digits=1)) \n" *
                "MLT $(round(longitude_to_MLT(ϕ_eq, degrees=false), digits=1)) \n",
    space     = :relative,
    align     = (:left, :top),
    fontsize  = 11
)

axislegend(ax, position=:rb)

#save("src/WPI/test_scripts/resonance/resonance_lat_vs_energy_log.png", fig)
