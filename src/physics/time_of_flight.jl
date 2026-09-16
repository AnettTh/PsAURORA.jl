using AURORA
using AURORA; mₑ, eV_in_J, c₀
using QuadGK
using LinearAlgebra

#===================================Saito-Miyoshi method===================================#

# IDEA: If these are slow, would it help to make them the same function but with multiple dispatch?
# NOTE: Verify that returning the absolute value is numerically allowed
"""
    t_whistler_transit(θ_resonance, v_g, R0)

Calculate transit time for whistler waves.

For some given group velocity of parallel propagating whistler waves, `v_g`, the time of
travel from the equatorial plane (`θ = π/2°`) to the resonance region `θ_resonance`, where θ
is the colatitude, is calculated.

# Arguments

- `θ_resonance`: Colatitude of the resonance region [rad].
- `v_g`: Group velocity of whistler waves [m/s].
- `R0`: Radius at the equatorial plane [m].

"""
function t_whistler_transit(θ_resonance, v_g, R0)

    function t(θ)
        return (R0 * sqrt(1 + 3*(cos(θ))^2) * sin(θ)) / v_g
    end

    t_w = quadgk(t, π/2, θ_resonance)

    return abs(t_w[1])
end


# TODO: Docstring
function t_electron_transit(θ_resonance, v_parallel_lc, R0)

    function t(θ)
        return (R0 * sqrt(1 + 3*(cos(θ))^2) * sin(θ)) / abs(v_parallel_lc)
    end

    t_e = quadgk(t, θ_resonance, π/2)

    return abs(t_e[1])
end


# TODO: Finish docstring
"""
    t_WPI_electron_precipitation(θ_resonance, v_g, v_parallel_lc, R0, t_l, E_eV, α_lc)

Calculate the total time-of-flight for a
Using the method of Saito-Miyoshi, as described in Saito 2012.
"""
function t_WPI_electron_precipitation(θ_resonance, v_g, v_parallel_lc, R0, t_l, E_eV, α_lc)

    t_w = t_whistler_transit(θ_resonance, v_g, R0)
    t_e = t_electron_transit(θ_resonance, v_parallel_lc, R0)
    t_b = quarter_bounceperiod(R0/RE, E_eV, mₑ, α_lc)

    t_precipitation = t_w + t_e + t_b + t_l

    return t_precipitation
end



#====================================Chen/Hsieh method====================================#

function chorus_angular_frequency_timer(ω)

end


function t_WPI_adiabatic()

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
