using AURORA

# IDEA: Make this into 'AbstractPropagation'
function time_of_flight(
    E,
    μ,
    z_distance;
    propagation=:simple,
    magnetic_field=nothing,
    r0=nothing,
    z_end=nothing
)

    if propagation == :simple
        println("Using simple propagation")
        return z_distance / (abs(μ) * v_of_E(E))

    elseif propagation == :fieldline

        # Construct local coordinate system based on magnetic field in initial position
        B = magnetic_field(r0...)
        b, e1, e2 = magnetic_basis(B)

        # Construct initial velocity with available information
        # TODO: Look into α_frac here, how do i know it is correct with regards to the rest?
        v0 = get_v0(dipole_field, r0, E, mₑ)

        result = boris_mover_TOF(magnetic_field, r0, v0, z_end)

        println("Using field-line propagation")
        return result.tof

    else
        throw(
            ArgumentError(
                "Non-valid propagation-mode chosen, use `:simple` or `:fieldline`"
                )
            )
    end
end
