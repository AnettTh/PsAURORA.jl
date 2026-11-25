using AURORA
using Dates

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
density_model = DentonDensity(fieldline)

## Define variables and make the distribution
N = electron_density(density_model, fieldline.position_eq)

T_parallel = 20e3
A = 0.0
Δ = 0.0
r = 2.0
M = 20
α_lc = losscone_angle(fieldline)

a_par = sqrt(2 * T_parallel * eV_in_J / mₑ)
v_par_grid  = range(-4.0 * a_par, 4.0 * a_par, length=300)
v_perp_grid = range(-4.0 * a_par, 4.0 * a_par, length=300)

F = [subtracted_bimaxwellian(vpa, vpe, N, T_parallel, A, α_lc, Δ, r, M)
    for vpa in v_par_grid, vpe in v_perp_grid]

## Frequencies to test
Ω_eq = minimum(plasma.Ω_e)
ω_frac = [0.38]
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
        (E_eV < E_min || E_eV > E_max) && continue

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
    xlabel = "v∥/c₀",
    ylabel = "v⊥/c₀",
    title  = "Resonant particles in subtracted bi-Maxwellian for ω=$(ω_frac[1])Ωₑ\n"*
    "T∥=$(round(T_parallel*1e-3, digits=1)) keV, A=$(round(A, digits=1)) keV, Δ=$Δ, r=$r",
    aspect=DataAspect()
)

## Meaningful range, taken only from the resonant cells
F_log  = log10.(max.(F, maximum(F) * 1e-6))
shown  = valid .& mask                              # tested and resonant
F_plot = ifelse.(shown, F_log, NaN)

any(shown) || error("No resonant cells, so there is no range to show")
clims = -19, maximum(F_plot[shown]) #extrema(F_plot[shown])

x = collect(v_par_grid) ./ c₀
y = collect(v_perp_grid) ./ c₀

hm = heatmap!(ax, x, y, F_plot';
    colormap   = :turbo,
    colorrange = clims,
    nan_color  = :white,
)
Colorbar(fig[1, 2], hm; label="log(Phase space density [m⁻⁶s³])", tellheight=false)
ylims!(ax, 0.0, 0.4)
xlims!(ax, -0.4, 0.4)
#save("src/WPI/test_scripts/distributions/figures/11232026ResonantDistribution$(ω_frac[1])Ωₑ.png", fig)
#save("src/WPI/test_scripts/distributions/figures/11232026ResonantDistribution$(ω_frac[1])Ωₑ_WithTreshold.png", fig)
