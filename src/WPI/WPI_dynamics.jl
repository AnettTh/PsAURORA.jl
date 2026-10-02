using AURORA
using AURORA; mₑ, eV_in_J, c₀
using QuadGK
using LinearAlgebra
using Roots

# IDEA: Make road-map of the WPI-file, as it is a bit messy at this point
"""
    dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)

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
- `ω_pe`: Electron plasma frequencie [rad/s].
- `Ω_e`: Electron cyclotron frequencie [rad/s].

# Returns

- The wave number that satisfies the dispersion relation.

# Throws

- `ArgumentError`: If the WNA is outside of ±1 (not in radians).
- `ArgumentError`: If the frequency is outside the LBC range.
"""
function dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)

    abs(θ) > 1 && throw(ArgumentError("Are you sure you are using radians for the WNA?"))

    any(ω .> abs.(0.5 .* Ω_e)) && throw(ArgumentError(
        "You are not in the LBC-range, sure this is right?"
        ))

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

# NOTE: Might need a WNA here?
"""
    group_velocity_whistler_wave(ω, Ω_e, ω_pe)

Calculates the group velocity `v_g` for the whistler mode chorus wave as a function of
chorus angular frequency `ω`.

This is valid assuming chorus frequencies `ω` ≫ ion gyro-frequencies (Chen 2020).

# Arguments

- `ω`: Chorus angular frequency, often a linearly rising tone [rad/s].
- `Ω_e`: Electron gyrofrequency at latitude λ [rad/s].
- `ω_pe`: Electron plasma frequency at latitude λ [rad/s].

# Returns

- The group velocity for the set of parameters chosen [m/s].

# Throws

- `ArgumentError`: If the frequencies does not match that of the whistler-branch.
"""
function group_velocity_whistler_wave(ω, Ω_e, ω_pe)

    any(@. ω > Ω_e) && throw(ArgumentError("Not on whistler branch, check your frequencies!"))

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

        ω_pe = field_dependent ?
            ωpe_at_λ(plasma.ne_model, particle.L, λ) : #, particle.ϕ, particle.magnetic_field) :
            plasma.ω_pe[1]

        v_g = group_velocity_whistler_wave(ω, Ω_e, ω_pe)
        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        return ds_dλ / v_g
    end

    t_w, _ = quadgk(f, 0.0, λ_resonance)
    return abs.(t_w)
end


"""
    particle_transit(particle, λ_resonance; z_ionosphere::Float64=600e3)

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
            # TODO: Also here hardcoded solution to fix
            if particle.magnetic_field == dipole_field
                B_λ = particle.magnetic_field(particle.L, λ)
            elseif particle.magnetic_field == tsyganenko_field
                B_λ = tsyganenko_field_spherical(particle.L, λ, particle.ϕ)
            else
                throw(ArgumentError("Hardcoded temporary solution, look into it!"))
            end

            α_λ = pitch_angle_at_λ(particle.α_eq, particle.B_eq, norm(B_λ))
            isnothing(α_λ) && return 0.0
            vz = particle.v * cos(α_λ)
            iszero(vz) && return 0.0
        else
            vz = abs(particle.v * cos(particle.α_lc))
        end

        # Removing NaN-values to not break the integration
        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        result = ds_dλ / vz
        isnan(result) && return 0.0
        return result
    end

    t_e, _ = quadgk(f, λ_resonance, particle.λ_ionosphere)
    return abs.(t_e)
end


# NOTE: This is not a realistic model, but used as a temporary solution. Look into default value of t!
"""
    wave_chirp(ω::AbstractVector; t::Real=0.2)

Compute the chorus wave launch time profile t₀(ω) at the magnetic equator.

Returns the time at which a wave of angular frequency ω is launched from the equatorial
source region, assuming a linear frequency chirp (rising tone), as modeled in Chen 2020. The
default parameter of rise-time in 0.2 s correspond to Chen 2020, and the start- and end
frequencies is the start- and end point of the frequency range.

# Arguments

- `ω`: Angular frequency at which to evaluate t₀ [rad/s].

# Keyword Arguments

- `t`: Duration of the chorus element [s], default is 0.2 s.
"""
function wave_chirp(ω::AbstractVector; t::Real=0.2)
    chirp_rate = ω[1] - ω[end] / t
    return (ω .- ω[1]) ./ chirp_rate
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

# Returns

- The latitude of resonance for the given frequency, particle and plasma [rad].
"""
function resonance_latitude(ω, particle, plasma; θ::Float64=0.0, n::Int=1)

    # Define the resonance condition (equation that should equal zero)
    function resonance_condition(λ; n=n, θ=θ)
        Ω_e = Ωe_at_λ(λ, particle.L, particle.magnetic_field)
        ω_pe = ωpe_at_λ(plasma.ne_model, particle.L, λ) #, particle.ϕ, particle.magnetic_field)

        k = dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)
        k_parallel = k * cos(θ)

        # TODO: remove hard-coded solution
        if particle.magnetic_field == dipole_field
            B = norm(particle.magnetic_field(particle.L, λ))
        elseif particle.magnetic_field == tsyganenko_field
            B = norm(tsyganenko_field_spherical(particle.L, λ, deg2rad(120)))
        else
            throw("Hardcoded temporary solution, check it!")
        end

        α = pitch_angle_at_λ(particle.α_eq, particle.B_eq, B)

        v_parallel = particle.v * cos(α)

        return ω - k_parallel * v_parallel + (n*Ω_e / particle.γ)
    end

    # Find which λ causes the resonance condition-function to change sign
    λ_grid = plasma.λ
    f_possible = resonance_condition.(λ_grid, n=n)

    idx = findfirst(i -> f_possible[i] * f_possible[i+1] < 0, 1:length(f_possible)-1)

    # Return nothing if there is no resonance
    #isnothing(idx) && return nothing
    isnothing(idx) && return NaN

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
    ω_grid::AbstractVector,
    particle::ParticleState,
    plasma::PlasmaState;
    wave_launch_time=nothing,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1
)

    WPI_TOF(ω::Real, particle, plasma; kwargs...

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
  from the source region (equator). The default is a simple chirp, but can be any function
  that takes frequencies `ω`.
- `field_dependent`: Decides which model to use, default is the field-dependent one.
- `θ`: Wave-normal angle of the wave, default is `1.0` (field-aligned).
- `n`: Harmonic number, default is `1`.

# Returns

- The time of flight for the particle in the plasma-state under the influence of the wave
  [s].
"""
function WPI_TOF(
    ω_grid::AbstractVector,
    particle::ParticleState,
    plasma::PlasmaState;
    wave_launch_time=wave_chirp,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1
)

    # Calculate wave-launch time based on given model
    # NOTE: Need to incorporate kwargs here in some way
    t_l = wave_launch_time(ω_grid)

    tof = zeros(length(ω_grid))

    # Loop over all frequencies as they have different resonance regions
    for (i, ω) in enumerate(ω_grid)
        λ_res = resonance_latitude(ω, particle, plasma; θ=θ, n=n)

        # For conbinations that don't resonate
        if isnothing(λ_res) || isnan(λ_res)
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
        # IDEA: Add boris-mover in magnetic field here!
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
