using AURORA


"""
    get_v0_from_Eμ(
    magnetic_field,
    r0,
    E_eV,
    μ0;
    ϕ=0.0,
    towards_equator::Bool=true,
    relativistic::Bool=false
)

Calculate the electron velocity vector given initial values.

Using the given magnetic field model, initial position, particle energy and pitch-angle to
find the appropriate initial velocity for the electron. To choose the direction of the
resulting velocity (i.e. away from or towards the Earth), specify `flip`, `false` being
towards the Earth, `true` being away from. If needed, it is also possible to control the
phase of the gyration.

# Arguments

- `magnetic_field`: Function describing the magnetic field, taking three positional values
  as the argument.
- `r_0`: The initial position of the particle [m].
- `E_eV`: Energy of the particle [eV].
- `μ0`: The cosine of the pitch-angle `α` of the particle, i.e. some value between -1 and 1.

# Keyword Arguments

- `ϕ`: The phase of the gyration, default is 0.0.
- `towards_equator`: Decides if the valocity is parallel or anti-parallel with the magnetic
  field, parallel being towards the equator IF in the northern hemisphere. The default is
  `true`.
- `relativistic`: Option to correct for relativistic effects, default is `false`.
"""
function get_v0_from_Eμ(
    magnetic_field,
    r0,
    E_eV,
    μ0;
    ϕ=0.0,
    towards_equator::Bool=true,
    relativistic::Bool=false
)

    -1 ≤ μ0 ≤ 1 || throw(ArgumentError("μ0 must be between -1 and 1"))

    v_mag = velocity_from_kinetic_energy(E_eV, mₑ; relativistic=relativistic)

    # Construct the magnetic basis
    B = magnetic_field(r0...)

    b, e1, e2, _ = magnetic_basis(B)

    # Adjust direction based on which hemisphere `r0` is in, meaning that if we are in the
    # northern hemisphere, towards equator means μ0 > 0
    if towards_equator
        μ0_corrected = r0[3] > 0 ? -abs(μ0) : abs(μ0)
    else
        μ0_corrected = r0[3] > 0 ? abs(μ0) : -abs(μ0)
    end

    # Determine parallel and perpendicular speeds
    v_parallel = μ0_corrected * v_mag
    v_perp = sqrt(1 - μ0_corrected^2) * v_mag

    # Initial velocity
    v0 = v_parallel .* b .+
         v_perp .* (cos(ϕ) .* e1 .+ sin(ϕ) .* e2)

    return Tuple(v0)
end


#===================================Run preferred method===================================#
# IDEA: Rewrite this to also fit with the WPI-versions?
"""
    time_of_flight(
    E_eV,
    μ0,
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
- `μ0`: Cosine of particle pitch-angle.
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
    μ0,
    z_distance;
    propagation=:simple,
    magnetic_field=nothing,
    r0=nothing,
    r_source=nothing
)

    E_eV > 0 || throw(ArgumentError("E_eV must be nonzero and positive"))

    if propagation == :simple
        return z_distance / (abs(μ0) * v_of_E(E_eV))

    elseif propagation == :fieldline

        if any(isnothing, (magnetic_field, r0, r_source))
            throw(ArgumentError(
                "For field-line propagation, `magnetic_field`, `r0` and `r_source` are " *
                "required. Check the keyword arguments for the function."
                )
            )
        end

        # Construct initial velocity with available information
        v0 = get_v0_from_Eμ0(magnetic_field, r0, E_eV, μ0; towards_equator=true)

        result = boris_mover_TOF(magnetic_field, r0, v0, r_source)


        isnothing(result) && throw(
            ArgumentError("Particle did not precipitate, i.e. no TOF available.")
        )

        return result.tof

    # NOTE: This is a placeholder! Remove or change.
    elseif propagation == :WPI_field_independent
        tof = 1.0

    elseif propagation == :WPI_field_dependent
        tof = 1.0

    else
        throw(
            ArgumentError(
                "Non-valid propagation-mode chosen, use `:simple`, `:fieldline`,
                `:WPI_field_independent` or `:WPI_field_dependent`."
                )
            )
    end
end
