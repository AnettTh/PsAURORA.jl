using AURORA
using AURORA; mₑ, eV_in_J, c₀
using QuadGK
using LinearAlgebra
using Roots


"""
    dispersion_relation_whistler_branch(ω, θ; ω_pe=37.9e3 * 2π, Ω_e=9.48e3 * 2π)

Calculate the wavenumber the whistler-branch of the Appleton-Hartree equation for oblique
waves (Hsieh 2022, equation ).

The function uses a defined whistler angular frequency `ω` and a wave-normal angle (WNA) `θ`
to solve the dispersion relation. The electron plasma frequency `ω_pe` and the electron
cyclotron frequency `Ω_e` are functions of position, while `ω` is a function of time,
meaning that the function supports ranges either as a funciton of position or time, but not
both.

# Arguments

- `ω`: Wave angular frequency for the whistler-mode chorus wave [rad/s].
- `θ`: Wave normal angle of the propagating wave, non-zero value indicates oblique wave
  [rad].

# Keyword Arguments

- `ω_pe`: Electron plasma frequencie, default is 37.9 kHz (Hsieh 2022) [Hz].
- `Ω_e`: Electron cyclotron frequencie, default is 9.48 kHz (Hsieh 2022) [Hz].

# Returns

- The wave number that satisfies the dispersion relation.

# Throws

- `ArgumentError`: If the WNA is outside of ±1 (not in radians).
"""
function dispersion_relation_whistler_branch(ω, θ; ω_pe=37.9e3*2π, Ω_e=9.48e3*2π)

    abs(θ) > 1 && throw(ArgumentError("Are you sure you are using radians for the WNA?"))

    ω > abs(0.5 * Ω_e) && throw(ArgumentError(
        "You are not in the LBC-range, sure this is right?"
        ))

    # TODO: Add cold-plasma check, I think that involves Ω_e and ω_pe??

    # NOTE: Can remove the throws for array-check when the rest of the code is safe
    if ω isa AbstractArray && Ω_e isa AbstractArray
        length(ω) == length(Ω_e) || throw(DimensionMismatch(
            "ω and Ω_e must have the same length when both are arrays — " *
            "they represent different physical dimensions (frequency vs latitude)"
        ))
    end
    if ω isa AbstractArray && ω_pe isa AbstractArray
        length(ω) == length(ω_pe) || throw(DimensionMismatch(
            "ω and ω_pe must have the same length when both are arrays"
        ))
    end

    X = @. ω_pe^2 / ω^2
    Y = @. Ω_e / ω

    a = @. 2*(1-X)
    sin2 = sin(θ)^2

    num = @. X*a
    denum = @. a - Y^2 * sin2 + Y*sqrt(Y^2 * sin2^2 + a^2 * cos(θ)^2)

    c2k2ω2 = @. 1 - (num/denum)

    k = @. sqrt(c2k2ω2) * ω / c₀

    return k
end


"""
    group_velocity_whistler_wave(
    ω::AbstractArray;
    Ω_e::Float64=9.48e3 * 2π,
    ω_pe::Float64=37.9e3 * 2π)

Calculates the group velocity `v_g` for the whistler mode chorus wave as a function of
chorus angular frequency `ω`.

This is valid assuming chorus frequencies `ω` ≫ ion gyro-frequencies (Chen 2020).

# Arguments

- `ω`: Chorus angular frequency, often a linearly rising tone [rad/s].

# Keyword Arguments

- `Ω_e`: Electron gyrofrequency at a specific location λ [rad/s], default is 9.48 kHz (from
  Hsieh et al. 2022).
- `ω_pe`: Electron plasma frequency at a specific location λ [rad/s], default is
  `4×Ω_e`=37.9 kHz (from Hsieh et al. 2022).

# Returns

- The group velocity for the set of parameters chosen.

# Throws

- `ArgumentError`: If the frequencies does not match that of the whistler-branch.
- `ArgumentError`: If the plasma parameters are outside of the expected domain, might not be
  non-realistic, but for now it indicates an error.
"""
function group_velocity_whistler_wave(
    ω::Union{Real, AbstractArray};
    Ω_e::Union{Float64, AbstractVector}=9.48e3*2π,
    ω_pe::Union{Float64, AbstractVector}=37.9e3*2π)

    # TODO: add cold-plasma check?
    any(@. ω > Ω_e) && throw(ArgumentError("Not on whistler branch, check your frequencies!"))
    #any(ω_pe .< Ω_e) && throw(ArgumentError("Plasma parameters not valid, reality-check needed!"))

    a = @. (2 * c₀) / (ω_pe / Ω_e)
    b = @. (1 - (ω / Ω_e))^(3/2)
    c = @. (ω / Ω_e)^(1/2)

    return @. a*b*c
end


"""
    wave_transit(
    ω,
    λ_resonance,
    particle::ParticleState,
    plasma::PlasmaState;
    field_dependent::Bool=true)

Calculates the transit time of the whistler wave from the equator to the resonance latitude.

By either using a field-independent model (Miyoshi and Saito, 2012) or a field-dependent
model (Hsieh et al. 2022), where the latter is default. The function calculates the time it
takes for the wave specified by the frequency `ω`, the particle- and plasma state, traveling '
from the source region (equator) to the resonance latitude `λ_resonance`.

# Arguments

- `ω`: Wave angular frequency of the whistler wave [rad/s].
- `λ_resonance`: The latitude of the resonance region (0.0 is equator) [rad].
- `particle`: Structure holding the particle state.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `field_dependent`: Decides which model to use, default is the field-dependent one.

# Returns

- Wave transit time for each frequency ω [s].
"""
function wave_transit(
    ω,
    λ_resonance,
    particle::ParticleState,
    plasma::PlasmaState;
    field_dependent::Bool=true)

    R0 = particle.L * RE

    function f(λ)
        # If field-dependent, choose based on position, else use equatorial Ω_e
        Ω_e = field_dependent ?
            Ωe_at_λ(λ, particle.L, particle.magnetic_field) :
            plasma.Ω_e[1]

        # NOTE: This currently assumes nₑ constant along the field-line. Look for model?
        v_g = group_velocity_whistler_wave(ω; Ω_e=Ω_e, ω_pe=plasma.ω_pe[1])
        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        return ds_dλ / v_g
    end

    t_w, _ = quadgk(f, 0.0, λ_resonance)
    return abs.(t_w)
end


"""
    particle_transit(
    particle,
    λ_resonance;
    z_ionosphere::Float64=600e3
    )

Calculates the transit time of the particle from the resonance latitude to the ionosphere.

By either using a field-independent model (Miyoshi and Saito, 2012) or a field-dependent
model (Hsieh et al. 2022), where the latter is default. The function calculates the time it
takes for the particle in the defined state to travel from the resonance region
`λ_resonance` to the defined ionospheric latitude, `λ_ionosphere`, as defined in the
particle state (600 km altitude).

# Arguments

- `particle`: Structure holding the particle state.
- `λ_resonance`: The latitude of the resonance region (0.0 is equator) [rad].

# Keyword Arguments

- `field_dependent`: Decides which model to use, default is the field-dependent one.

# Returns

- Particle transit time for each frequency ω [s].
"""
function particle_transit(particle, λ_resonance; field_dependent::Bool=true)

    R0 = particle.L * RE

    function f(λ)
        if field_dependent
            B_λ = particle.magnetic_field(particle.L, λ)
            α_λ = pitch_angle_at_λ(particle.α_eq, particle.B_eq, norm(B_λ))
            isnothing(α_λ) && return 0.0
            vz = particle.v * cos(α_λ)
            iszero(vz) && return 0.0
        else
            vz = abs(particle.v * cos(particle.α_lc))       # NOTE: becomes relativistic from particle.v, if needed
        end

        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        return ds_dλ / vz
    end

    # NOTE: will not run if the particle is outside the loss-cone, as the limits here is invalid, but that might be fine? Add check for this?
    t_e, _ = quadgk(f, λ_resonance, particle.λ_ionosphere)


    return abs.(t_e)
end


function wave_chirp(ω; ω0=2π*600, ω1=2π*1350, t=0.2)
    chirp_rate = (ω0 - ω1) / t
    return (ω .- ω0) ./ chirp_rate
end


"""
    resonance_latitude(ω, particle, plasma; θ::Float64=0.0, n::Int=1)

Calculate the resonance latitude for the defined particle with a whistler mode wave.

Depending on the frequency of the wave, the state of the particle and the state of the
plasma, the resonance latitude in a magnetic field is found. `n` denotes the harmonic, `1`
being the default value.

# Arguments

- `ω`: The frequency of the wave [rad/s].
- `particle`: Structure holding the state of the particle.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `θ`: WNA, default is nothing (ducted).
- `n`: Resonance number, default is `1` (cyclotron resonance).
- `relativistic`: Option to correct for relativistic energies, default is `false`.

# Returns

- The latitude of resonance for the given frequency, particle and plasma [rad].
"""
function resonance_latitude(ω, particle, plasma; θ::Float64=0.0, n::Int=1)

    # TODO: Change this to work also for ω_pe as a function of position
    ω_pe = plasma.ω_pe[1]  # constant

    # Define the resonance condition (equation that should equal zero)
    function resonance_condition(λ; n=n, θ=θ)
        # NOTE: This already exist in plasma???
        Ω_e = Ωe_at_λ(λ, particle.L, particle.magnetic_field)

        k = dispersion_relation_whistler_branch(ω, θ; ω_pe=ω_pe, Ω_e=Ω_e)
        k_parallel = k * cos(θ)

        B = norm(particle.magnetic_field(particle.L, λ))
        α = pitch_angle_at_λ(particle.α_eq, particle.B_eq, B)

        v_parallel = particle.v * cos(α)

        return ω - k_parallel * v_parallel + (n*Ω_e / particle.γ)
    end

    # Find which λ causes the resonance condition-function to change sign
    λ_grid = plasma.λ
    f_possible = resonance_condition.(λ_grid, n=n)

    idx = findfirst(i -> f_possible[i] * f_possible[i+1] < 0, 1:length(f_possible)-1)

    # Return nothing if there is no resonance
    isnothing(idx) && return nothing

    # Figure out the latitude where the funciton changed sign
    λ_resonance = find_zero(resonance_condition, (λ_grid[idx], λ_grid[idx+1]))

    return λ_resonance
end

function resonance_latitude(ω_grid::AbstractVector, particle, plasma; θ::Float64=0.0, n::Int=1)
    return [let r = resonance_latitude(ω, particle, plasma; θ=θ, n=n)
                isnothing(r) ? NaN : r
            end for ω in ω_grid]
end


"""
    WPI_TOF(
    ω_grid,
    particle,
    plasma;
    wave_launch_time=nothing,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1
)

Calculates the time-of-flight of a particle resonating with a whistler mode chorus wave,
then precipitating into the ionosphere.

The funciton uses a defined frequency-grid for the wave, particle state and plasma state to
calculate the resonance latitude. The wave launch time can be given as a function in the
keyword arguments, else, it is assumed instant at all frequencies.

# Arguments

- `ω_grid`: The frequency grid of the wave [rad/s].
- `particle`: Structure holding the state of the particle.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `wave_launch_time`: The time at which the frequencies in the frequency grid is launched
  from the source region (equator). If everything is launched at once, this is set to
  `nothing` (default), or it can be a defined function that takes frequencies `ω`.
- `field_dependent`: Decides which model to use, default is the field-dependent one.
- `θ`: Wave-normal angle of the wave, default is `1.0` (field-aligned).
- `n`: Harmonic number, default is `1`.

# Returns

- The time of flight for the particle in the plasma-state under the influence of the wave
  [s].
"""
function WPI_TOF(
    ω_grid,
    particle,
    plasma;
    wave_launch_time=nothing,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1
)

    # Is either a scalar (Saito-Miyoshi) or a functon of frequency (Chen)
    t_l = isnothing(wave_launch_time) ? 0.0 : wave_launch_time(ω_grid)

    tof = zeros(length(ω_grid))

    # Loop over all frequencies as they have different resonance regions
    for (i, ω) in enumerate(ω_grid)
        λ_res = resonance_latitude(ω, particle, plasma; θ=θ, n=n)

        # For conbinations that don't resonate
        if isnothing(λ_res)
            tof[i] = NaN
            continue
        end

        # Calculate the components for the given ω and sum them
        t_w = wave_transit(
            ω,
            λ_res,
            particle,
            plasma;
            field_dependent=field_dependent
        )
        t_e = particle_transit(
            particle,
            λ_res;
            field_dependent=field_dependent
        )

        tof[i] = t_w + t_e + t_l[i]

    end
    return tof
end

function WPI_TOF(ω::Real, particle, plasma; kwargs...)
    return WPI_TOF([ω], particle, plasma; kwargs...)[1]
end
