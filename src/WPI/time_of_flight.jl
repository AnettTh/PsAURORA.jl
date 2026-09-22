using AURORA
using AURORA; mₑ, eV_in_J, c₀
using QuadGK
using LinearAlgebra
using Roots


"""
    group_velocity_whistler_wave(
    ω::AbstractArray;
    Ω_e::Float64=9.48e3,
    ω_pe::Float64=37.9e3)

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
    Ω_e::Union{Float64, AbstractVector}=9.48e3,
    ω_pe::Union{Float64, AbstractVector}=37.9e3)

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

- `ω`: Wave angular frequency of the whistler wave [rad/s]. # TODO: check if it is in Hz somewhere, and update docs!!!
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
            vz = abs(particle.v * cos(particle.α_lc))
        end

        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        return ds_dλ / vz
    end

    # NOTE: will not run if the particle is outside the loss-cone, as the limits here is invalid, but that might be fine?
    t_e, _ = quadgk(f, λ_resonance, particle.λ_ionosphere)


    return abs.(t_e)
end


# TODO: Add model for reaching the desired frequency (linear rising tone?)
function wave_launch_time(ω; ω0=2π*600, ω1=2π*1350, t=0.2)
    chirp_rate = (ω0 - ω1) / t
    return (ω .- ω0) ./ chirp_rate
end


# To find the point of resonance (eq. 6)
# TODO: Add γ and θ as kwarg???
# TODO: Replace this with the full dispersion-relation
"""
    parallel_wavenumber(ω_pe, Ω_e, ω)

Calculate the parallel wave number that satisfies the dispersion relation.
"""
function parallel_wavenumber(ω_pe, Ω_e, ω)

    frac = @. ω_pe^2 / (ω*(abs(Ω_e) - ω))
    k_parallel = @. (ω / c₀) * sqrt(1 + frac)

    return k_parallel
end



# TODO: Compute this once, and remove in-function calculations in the other functions, as there should be a yes/no option for relativistic particles for comparison
function γ(E_eV)
    return 1
end


# TODO: Add γ and θ as optional parameters!
"""
    resonance_latitude(ω, particle, plasma; γ=nothing, θ=nothing, n::Int=1)

Calculate the resonance latitude for the defined particle with a whistler mode wave.

Depending on the frequency of the wave, the state of the particle and the state of the
plasma, the resonance latitude in a magnetic field is found. `n` denotes the harmonic, `1`
being the default value.

# Arguments

- `ω`: The frequency of the wave [rad].
- `particle`: Structure holding the state of the particle.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `γ`
- `θ`
- `n`: Resonance number, default is `1` (cyclotron resonance).

# Returns

- The latitude of resonance for the given frequency, particle and plasma [rad].
"""
function resonance_latitude(ω, particle, plasma; γ=nothing, θ=nothing, n::Int=1)

    # TODO: Change this to work also for ω_pe as a function of position
    ω_pe = plasma.ω_pe[1]  # constant

    # Define the resonance condition (equation that should equal zero)
    function resonance_condition(λ; n=n)
        # NOTE: This already exist in plasma???
        Ω_e = Ωe_at_λ(λ, particle.L, particle.magnetic_field)

        k = parallel_wavenumber(ω_pe, Ω_e, ω)

        B = norm(particle.magnetic_field(particle.L, λ))
        α = pitch_angle_at_λ(particle.α_eq, particle.B_eq, B)

        v = particle.v * cos(α)

        return ω - k * v + n*Ω_e
    end

    # Find which λ causes the resonance condition-function to change sign
    λ_grid = plasma.λ
    f_possible = resonance_condition.(λ_grid)

    idx = findfirst(i -> f_possible[i] * f_possible[i+1] < 0, 1:length(f_possible)-1)

    # Return nothing if there is no resonance
    isnothing(idx) && return nothing

    # Figure out the latitude where the funciton changed sign
    λ_resonance = find_zero(resonance_condition, (λ_grid[idx], λ_grid[idx+1]))

    return λ_resonance
end


function resonance_latitude(ω_grid::AbstractVector, particle, plasma; γ=nothing, θ=nothing, n::Int=1)
    return [resonance_latitude(ω, particle, plasma; γ=γ, θ=θ, n=n) for ω in ω_grid]
end


# TODO: Finish docstring
"""
    WPI_TOF(
    ω,
    λ_resonance,
    particle,
    plasma;
    wave_launch_time=nothing,
    field_dependent::Bool=true,
    θ::Float64=1.0
)

Calculates the time-of-flight of a particle resonating with a whistler mode chorus wave,
then precipitating into the ionosphere.

The funciton uses a defined frequency-grid for the wave, particle state and plasma state to
calculate the resonance latitude. The wave launch time can be given as a function in the
keyword arguments, else, it is assumed instant at all frequencies.

# Arguments

- `ω_grid`: The frequency grid of the wave [rad/s].
- `λ_resonance`: The latitude of resonance [rad].
- `particle`: Structure holding the state of the particle.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `wave_launch_time`: The time at which the frequencies in the frequency grid is launched
  from the source region (equator). If everything is launched at once, this is set to
  `nothing` (default), or it can be a defined function that takes frequencies `ω`.
- `field_dependent`: Decides which model to use, default is the field-dependent one.
- `θ`: Wave-normal angle of the wave, default is `1.0` (field-aligned).

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
    θ::Float64=1.0
)

    # Is either a scalar (Saito-Miyoshi) or a functon of frequency (Chen)
    t_l = isnothing(wave_launch_time) ? 0.0 : wave_launch_time(ω_grid)

    tof = zeros(length(ω_grid))

    # Loop over all frequencies as they have different resonance regions
    for (i, ω) in enumerate(ω_grid)
        λ_res = resonance_latitude(ω, particle, plasma, θ=θ)

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


#===================================Run preferred method===================================#
# IDEA: Make this into 'AbstractPropagation'
# TODO: z_distance only used if simple, is there a better way? Change to kwarg with added errors?
# TODO: Add TOF including also whistler wave influence
"""
    time_of_flight(
    E_eV,
    μ,
    z_distance;
    propagation=:simple,
    magnetic_field=nothing,
    r0=nothing,
    r_source=nothing
)

Calculate time-of-flight (TOF) for a particle given its energy and pitch-angle.

If `propagation=:simple`, the TOF is approximated by the distance traveled, `z_distance`
[m], divided by the component of the velocity parallel to the magnetic field [m/s]. For
`propagation=:fieldline`, a boris-mover numerical scheme is used to simulate the time of
flight.

# Arguments

- `E_eV`: Particle energy (eV).
- `μ`: Cosine of particle pitch-angle.
- `z_distance`: Distance traveled before hitting the ionosphere [m]. ONLY USED FOR THE SIMPLE PROPAGATION MODE!!

# Keyword Arguments

- `propagation`: Propagation mode, either `:simple` or `:fieldline`. Default is `:simple`.
- `magnetic_field`: Magnetic field function `f(x, y, z)`. Required for `:fieldline`.
- `r0`: Initial position vector [m]. Required for `:fieldline`.
- `r_source`: Stopping radial distance from Earth's center [m]. Required for `:fieldline`.

# Returns

- The time of flight using the chosen propagation method.

# Throws

- `ArgumentError`: If `E_eV` is zero or negative.
- `ArgumentError`: If `propagation` is not `:simple` or `:fieldline`.
- `ArgumentError`: If `propagation=:fieldline` and any of `magnetic_field`, `r0`, or
  `r_source` are `nothing`.
- `ArgumentError`: If the particle does not precipitate and no time-of-flight is available.
"""
function time_of_flight(
    E_eV,
    μ,
    z_distance;
    propagation=:simple,
    magnetic_field=nothing,
    r0=nothing,
    r_source=nothing
)

    E_eV > 0 || throw(ArgumentError("E_eV must be nonzero and positive"))

    if propagation == :simple
        return z_distance / (abs(μ) * v_of_E(E_eV))

    elseif propagation == :fieldline

        if any(isnothing, (magnetic_field, r0, r_source))
            throw(ArgumentError(
                "For field-line propagation, `magnetic_field`, `r0` and `r_source` are " *
                "required. Check the keyword arguments for the function."
                )
            )
        end

        # Construct initial velocity with available information
        v0 = get_v0_from_Eμ(magnetic_field, r0, E_eV, μ; towards_equator=true)

        result = boris_mover_TOF(magnetic_field, r0, v0, r_source)


        isnothing(result) && throw(
            ArgumentError("Particle did not precipitate, i.e. no TOF available.")
        )

        return result.tof

    else
        throw(
            ArgumentError(
                "Non-valid propagation-mode chosen, use `:simple` or `:fieldline`"
                )
            )
    end
end
