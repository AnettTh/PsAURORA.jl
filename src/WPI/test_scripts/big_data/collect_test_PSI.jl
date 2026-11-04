using AURORA
using JLD2
using ProgressMeter

## Magnetic field and fixed settings
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
use_boris = true

## Grids (dimension order of TOF: E, μ, ω, L, Ψ0)
L_grid = [4.0, 6.0]
E_grid = range(1e3, 40e3, length=10)
α_grid = [deg2rad(0.5), deg2rad(2)]
μ_grid = -cos.(α_grid)
ω_frac = range(0.1, 0.5, length=10)
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
Ψ0_grid = collect(0:π/4:7π/4)

## Wave amplitude is constant for all waves
amplitude = ConstantAmplitude(100e-12)

## Output
TOF = fill(NaN, length(E_grid), length(μ_grid), length(ω_frac), length(L_grid), length(Ψ0_grid))
Ω_e_eq = zeros(length(L_grid))
done = falses(length(L_grid), length(Ψ0_grid))
outfile = "src/WPI/data/02112026_Psi0test_MLT$(round(longitude_to_MLT(ϕ_eq; degrees=false), digits=1)).jld2"

function save_results(file, TOF, Ω_e_eq, done)
    jldsave(file;
        TOF, Ω_e_eq, done,
        E_grid, μ_grid, α_grid, ω_frac, L_grid, Ψ0_grid,
        ϕ_eq, Bw = amplitude.Bw, use_boris,
    )
end

## Compute
for (iL, L) in enumerate(L_grid)
    # Make initial position used to define particles
    r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))

    # Make the plasma
    plasma = PlasmaState(λ_grid, ϕ_eq, magnetic_field, L)

    # Make the wave
    Ω_e_eq[iL] = minimum(plasma.Ω_e)
    ωs = ω_frac .* Ω_e_eq[iL]
    modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

    #
    pairs = CartesianIndices((length(E_grid), length(μ_grid)))

    # Particles depend on (E, μ, L) but not on Ψ0, so build them once per L
    particles = Matrix{Union{Nothing, ParticleState}}(nothing, size(pairs)...)
    Threads.@threads :greedy for k in eachindex(pairs)
        iE, iμ = Tuple(pairs[k])
        try
            particles[iE, iμ] = ParticleState(
                E_grid[iE],
                μ_grid[iμ],
                r0,
                magnetic_field,
                ϕ_eq;
                relativistic=true,
            )
        catch e
            @warn "ParticleState failed" L E=E_grid[iE] μ=μ_grid[iμ] exception=e
        end
    end

    # Running the simulation for each wave- and particle combination
    for (iΨ, Ψ0) in enumerate(Ψ0_grid)
        # Make wave with current Ψ0
        wave = WaveState(modulation, amplitude, plasma; Ψ0=Ψ0)

        # Track the progress
        progress = Progress(
            length(pairs);
            desc="L = $L, Ψ0 = $(round(rad2deg(Ψ0), digits=0))°: ",
        )
        n_failed = Threads.Atomic{Int}(0)

        # Run the TOF calculation for each of the energy- and pitch-angle pairs
        Threads.@threads :greedy for k in eachindex(pairs)
            iE, iμ = Tuple(pairs[k])
            p = particles[iE, iμ]
            if !isnothing(p)
                try
                    TOF[iE, iμ, :, iL, iΨ] = WPI_TOF(
                        ωs,
                        p,
                        plasma;
                        wave=wave,
                        use_boris=use_boris,
                        field_dependent=true,
                    )
                catch e
                    # Print only the first failure, count the rest
                    if Threads.atomic_add!(n_failed, 1) == 0
                        @error "WPI_TOF failed" L Ψ0 iE iμ exception=(e, catch_backtrace())
                    end
                end
            end
            next!(progress)
        end

        n_failed[] > 0 && @warn "L = $L, Ψ0 index $iΨ: $(n_failed[]) pairs failed"
        done[iL, iΨ] = true

        # A checkpoint after every (L, Ψ0)
        save_results(outfile, TOF, Ω_e_eq, done)
    end
end
