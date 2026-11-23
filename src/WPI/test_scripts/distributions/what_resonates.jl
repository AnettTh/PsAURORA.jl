using AURORA
using Dates
# TODO: Look into these variables

## Define fieldline, plasma and wave to find resonance
t = DateTime("2020-12-21T02:00:00")
magnetic_field = TsyganenkoMagneticField(t; pdyn=2.0, dst=-28.0, byimf=0.0, bzimf=2.0)
fieldline = StationFieldLine(magnetic_field, LTPReference(:sodankyla), t)
plasma = PlasmaState(fieldline, DentonDensity(fieldline))
wave = WaveState(
    ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2),
    ConstantAmplitude(100e-12),
    plasma,
    fieldline;
    Ψ0=0.0
)

## Define variables and make the distribution
T_parallel = 10e3   # Most likely not correct
N = 1e6             # Get this from fieldline
A = 1.0
Δ = 0.2
r = 2.0
M = 20
α_lc = deg2rad(1.1) # Get this from fieldline?

a_par = sqrt(2 * T_parallel * eV_in_J / mₑ)
v_par_grid  = range(-5.0 * a_par, 5.0 * a_par, length=300)
v_perp_grid = range(-5.0 * a_par, 5.0 * a_par, length=300)

F = [subtracted_bimaxwellian(vpa, vpe, N, T_parallel, A, α_lc, Δ, r, M)
    for vpa in v_par_grid, vpe in v_perp_grid]

## Frequencies to test
Ω_eq = minimum(plasma.Ω_e)
ω_frac = [0.5]
ωs = ω_frac .* Ω_eq

## Resonance mask with the same shape as F
function resonance_mask(v_par_grid, v_perp_grid, fieldline, plasma, wave, ωs; E_min=7e3, E_max=40e3)

    mask = zeros(Bool, length(v_par_grid), length(v_perp_grid))
    valid = zeros(Bool, size(mask))

    Threads.@threads :greedy for k in eachindex(mask)
        i, j = Tuple(CartesianIndices(mask)[k])
        v_par, v_perp = v_par_grid[i], v_perp_grid[j]
        v = hypot(v_par, v_perp)

        # Nothing to test, or beyond the model's range
        (v == 0 || v ≥ 0.99 * c₀) && continue

        # Kinetic energy that gives exactly this speed, relativistically
        γ = 1 / sqrt(1 - (v / c₀)^2)
        E_eV = (γ - 1) * mₑ * c₀^2 / abs(qₑ)

        # Skip energies above and below some threshold
        #(E_eV < E_min || E_eV > E_max) && continue

        # Equatorial pitch angle, with μ0 = -v∥/v (ParticleState uses v∥ = -v·μ0)
        μ0 = -v_par / v

        valid[i, j] = true
        p = ParticleState(E_eV, μ0, fieldline; relativistic=true)

        # Check for resonance
        for ω in ωs
            res = try
                point_of_resonance(ω, p, fieldline, plasma, wave)
            catch
                nothing
            end
            if !isnothing(res)
                mask[i, j] = true
                break
            end
        end
    end
    return mask, valid
end

mask, valid = resonance_mask(v_par_grid, v_perp_grid, fieldline, plasma, wave, ωs)

@show count(mask), count(valid)

## Make plot
fig = Figure(size=(700, 400))
ax  = Axis(fig[1, 1];
    xlabel = "v∥ / c",
    ylabel = "v⊥ / c",
    title  = "Resonant particles in subtracted bi-Maxwellian for ω=$(ω_frac[1])Ωₑ\n"*
    "T∥=$(round(T_parallel*1e-3, digits=1)) keV, A=$(round(A, digits=1)) keV, Δ=$Δ, r=$r",
    aspect=DataAspect()
)

## Meaningful range
F_log = log10.(max.(F, maximum(F) * 1e-6))
F_plot = ifelse.(valid, F_log, NaN)
clims = extrema(filter(isfinite, F_plot))

hm = heatmap!(ax,
    collect(v_par_grid) ./ c₀,
    collect(v_perp_grid) ./ c₀,
    F_plot';
    colormap = :turbo,
    colorrange = clims,
    nan_color = :white
)

Colorbar(fig[1, 2], hm; label="log(Phase space density [m⁻⁶s³])", tellheight=false)
contour!(ax, v_par_grid / c₀, v_perp_grid / c₀, Float64.(mask)'; levels=[0.5], color=:white, linewidth=2)
ylims!(ax, 0.0, 0.99)
save("src/WPI/test_scripts/distributions/figures/11232026ResonantDistribution$(ω_frac[1])Ωₑ.png", fig)
#save("src/WPI/test_scripts/distributions/figures/11232026ResonantDistribution$(ω_frac[1])Ωₑ_WithTreshold.png", fig)
