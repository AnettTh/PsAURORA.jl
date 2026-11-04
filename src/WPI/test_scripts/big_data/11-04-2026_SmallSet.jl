using AURORA
using JLD2
using ProgressMeter

## Magnetic field and fixed settings
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
use_boris = true
L = 6.5

## Grids with TOF(E, μ, ω, Ψ0) = (10, 8, 5, 4)
E_grid = range(25e3, 30e3, length=10)
α_grid = range(deg2rad(0.8), deg2rad(1.0), length=8)
μ_grid = -cos.(α_grid)
ω_frac = range(0.3, 0.35, length=5)
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
Ψ0_grid = collect(0:π/2:7π/4)

## Initial position and wave amplitude is constant for all runs
r0 = Cartesian(Spherical(L, 0.0, ϕ_eq))
amplitude = ConstantAmplitude(100e-12)


## Output
TOF = zeros(length(E_grid), length(μ_grid), length(ω_frac), length(Ψ0_grid))
outfile = "src/WPI/data/11-04-2026/SmallTestSet.jld2"

function save_results(file, TOF, Ω_e_eq)
    jldsave(
        file;
        TOF,
        Ω_e_eq,
        E_grid,
        α_grid,
        ω_frac,
        L,
        Ψ0_grid,
        ϕ_eq,
        Bw=amplitude.Bw,
        use_boris,
    )
end



## Make the plasma
plasma = PlasmaState(λ_grid, ϕ_eq, magnetic_field, L)

## Define the frequencies used
Ω_e_eq = minimum(plasma.Ω_e)
ωs = ω_frac .* Ω_e_eq

## Make the wave modulation for the full lower-band whistler mode chorus wave
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

## Make all combinations of energies and pitch-angles
index_pairs = CartesianIndices((length(E_grid), length(μ_grid)))

## Build the grid of particles for each combination of E and μ
particles = Matrix{Union{Nothing, ParticleState}}(nothing, size(index_pairs)...)
Threads.@threads :greedy for k in eachindex(index_pairs)
        iE, iμ = Tuple(index_pairs[k])
        particles[iE, iμ] = ParticleState(
            E_grid[iE],
            μ_grid[iμ],
            r0,
            magnetic_field,
            ϕ_eq;
            relativistic=true,
        )
end

## Running the simulation for each wave- and particle combination
for (iΨ, Ψ0) in enumerate(Ψ0_grid)
    # Make wave with current Ψ0
    wave = WaveState(modulation, amplitude, plasma; Ψ0=Ψ0)

    # Track the progress
    progress = Progress(length(index_pairs); desc="Ψ₀ = $(round(rad2deg(Ψ0), digits=0))°: ")

    # Run the TOF calculation for each of the energy- and pitch-angle pairs
    Threads.@threads :greedy for k in eachindex(index_pairs)
        iE, iμ = Tuple(index_pairs[k])
        p = particles[iE, iμ]
        TOF[iE, iμ, :, iΨ] = WPI_TOF(
            ωs,
            p,
            plasma;
            wave=wave,
            use_boris=use_boris,
            field_dependent=true,
        )
        next!(progress)
    end
end

save_results(outfile, TOF, Ω_e_eq)

## Inspect data
data = load("src/WPI/data/11-04-2026/SmallTestSet.jld2")

E_grid = data["E_grid"] ./ 1e3
α_grid = rad2deg.(data["α_grid"])
ω_frac = data["ω_frac"]
Ψ0_grid = data["Ψ0_grid"]
TOF = data["TOF"]

## Make one figure for each pitch-angle
for iα in eachindex(α_grid)
    # Make figure for energy
    fig = Figure()
    Label(fig[0, 1:3], "Time-of-flight for α = $(round(α_grid[iα], digits=3))°";
        fontsize  = 18,
        font      = :bold,
        tellwidth = false,
    )

    # One for each Ψ₀
    ax1 = Axis(fig[1,1], ylabel="E [keV]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[1]), digits=0))°")
    ax2 = Axis(fig[1,2], title="Ψ₀=$(round(rad2deg(Ψ0_grid[2]), digits=0))°")
    ax3 = Axis(fig[2,1], xlabel="ω/Ωₑ", ylabel="E [keV]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[3]), digits=0))°")
    ax4 = Axis(fig[2,2], xlabel="ω/Ωₑ", title="Ψ₀=$(round(rad2deg(Ψ0_grid[4]), digits=0))°")

    # Same scale for the colorbar
    vals = [v for k in 1:4 for v in TOF[:, 4, :, k] if !isnan(v) && v != 0.0]
    colorrange = (minimum(vals) - 0.1, maximum(vals))

    # Make heatmap
    hm1 = heatmap!(ax1, ω_frac, E_grid, TOF[:, iα, :, 1]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm2 = heatmap!(ax2, ω_frac, E_grid, TOF[:, iα, :, 2]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm3 = heatmap!(ax3, ω_frac, E_grid, TOF[:, iα, :, 3]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm4 = heatmap!(ax4, ω_frac, E_grid, TOF[:, iα, :, 4]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)

    # Add the colorbar
    Colorbar(fig[1:2, 3], hm1; label="TOF [s]")

    save("src/WPI/test_scripts/big_data/figures/energy/TOF_$(round(α_grid[iα], digits=3)).png", fig)
end

## Make one figure for each energy
for iE in eachindex(E_grid)
    ## Make figure for pitch-angle
    fig = Figure()
    Label(fig[0, 1:3], "Time-of-flight for E = $(round(E_grid[iE], digits=1)) keV";
        fontsize  = 18,
        font      = :bold,
        tellwidth = false,
    )
    # One for each Ψ₀
    ax1 = Axis(fig[1,1], ylabel="α [°]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[1]), digits=0))°")
    ax2 = Axis(fig[1,2], title="Ψ₀=$(round(rad2deg(Ψ0_grid[2]), digits=0))°")
    ax3 = Axis(fig[2,1], xlabel="ω/Ωₑ", ylabel="α [°]", title="Ψ₀=$(round(rad2deg(Ψ0_grid[3]), digits=0))°")
    ax4 = Axis(fig[2,2], xlabel="ω/Ωₑ", title="Ψ₀=$(round(rad2deg(Ψ0_grid[4]), digits=0))°")

    # Same scale for the colorbar
    vals = [v for k in 1:4 for v in TOF[5, :, :, k] if !isnan(v) && v != 0.0]
    colorrange = (minimum(vals) - 0.1, maximum(vals))

    # Make heatmap
    hm1 = heatmap!(ax1, ω_frac, α_grid, TOF[iE, :, :, 1]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm2 = heatmap!(ax2, ω_frac, α_grid, TOF[iE, :, :, 2]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm3 = heatmap!(ax3, ω_frac, α_grid, TOF[iE, :, :, 3]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)
    hm4 = heatmap!(ax4, ω_frac, α_grid, TOF[iE, :, :, 4]'; colorrange=colorrange, colormap=:turbo)#, nan_color=:transparent, lowclip=:transparent)

    # Add the colorbar
    Colorbar(fig[1:2, 3], hm1; label="TOF [s]")

    save("src/WPI/test_scripts/big_data/figures/pitch-angle/TOF_$(round(E_grid[iE], digits=1))keV.png", fig)
end
