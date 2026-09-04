using AURORA

# IDEA: Make this into 'AbstractPropagation'
# TODO: z_distance only used if simple, is there a better way?
# TODO: Look into renaming 'z_end', as this one propagates backwards
function time_of_flight(
    E_eV,
    μ,
    z_distance;
    propagation=:simple,
    magnetic_field=nothing,
    r0=nothing,
    z_end=nothing
)

    if propagation == :simple
        return z_distance / (abs(μ) * v_of_E(E_eV))

    elseif propagation == :fieldline

        # Construct initial velocity with available information
        v0 = get_v0_from_Eμ(magnetic_field, r0, E_eV, μ)

        result = boris_mover_TOF(magnetic_field, r0, v0, z_end)

        isnothing(result) && error(
            "Particle did not precipitate, i.e. no TOF available."
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
