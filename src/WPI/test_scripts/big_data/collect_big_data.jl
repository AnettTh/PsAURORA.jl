using AURORA
using JLD2
using ProgressMeter

## Magnetic field and fixed settings
ϕ_eq           = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)
use_boris      = true

## Grids (dimension order of TOF: E, μ, ω, L)
L_grid = [4.0, 5.0, 6.0, 7.0, 8.0]
E_grid = range(10e3, 40e3, length=10)
α_grid = range(deg2rad(0.01), deg2rad(10), length=10)
μ_grid = -cos.(α_grid)
ω_frac = range(0.1, 0.5, length=10)
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)

## Wave: fixed phase so runs are reproducible and comparable between L
amplitude = ConstantAmplitude(100e-12)
Ψ0 = 0.0

## Output
TOF    = fill(NaN, length(E_grid), length(μ_grid), length(ω_frac), length(L_grid))
Ω_e_eq = zeros(length(L_grid))
outfile = "src/WPI/data/02112026_TOF_4D_MLT$(round(longitude_to_MLT(ϕ_eq; degrees=false), digits=1)).jld2"

function save_results(file, TOF, Ω_e_eq, iL_done)
    jldsave(file;
        TOF, Ω_e_eq, iL_done,
        E_grid = collect(E_grid), μ_grid, α_grid = collect(α_grid),
        ω_frac = collect(ω_frac), L_grid, ϕ_eq,
        Bw = amplitude.Bw, Ψ0, use_boris,
    )
end

## Compute
for (iL, L) in enumerate(L_grid)
    r0     = Cartesian(Dipolar(L, 0.0, ϕ_eq))
    plasma = PlasmaState(λ_grid, ϕ_eq, magnetic_field, L)

    Ω_e_eq[iL] = minimum(plasma.Ω_e)
    ωs         = ω_frac .* Ω_e_eq[iL]

    modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
    wave       = WaveState(modulation, amplitude, plasma; Ψ0=Ψ0)

    pairs    = CartesianIndices((length(E_grid), length(μ_grid)))
    progress = Progress(length(pairs); desc="L = $L: ")
    n_failed = Threads.Atomic{Int}(0)

    Threads.@threads :greedy for k in eachindex(pairs)
        iE, iμ = Tuple(pairs[k])
        try
            p = ParticleState(
                E_grid[iE], μ_grid[iμ], r0, magnetic_field, ϕ_eq;
                relativistic=true,
            )
            TOF[iE, iμ, :, iL] = WPI_TOF(
                ωs, p, plasma;
                wave=wave, use_boris=use_boris, field_dependent=true,
            )
        catch
            Threads.atomic_add!(n_failed, 1)   # leaves NaN in this slot
        end
        next!(progress)
    end

    n_failed[] > 0 && @warn "L = $L: $(n_failed[]) (E, μ) pairs threw an error"
    save_results(outfile, TOF, Ω_e_eq, iL)    # checkpoint after every L
end
