"""
    pitchangle(v, B)

Calculate the pitch angle for a particle.

Given the velocity vector and the magnetic field vector of a particle, the function
calculates the pitch angle.

# Arguments

- `v`: Velocity vector of the test particle `[vx, vy, vz]` [m/s].
- `B`: Magnetic field vector `[Bx, By, Bz]` [T].

# Keyword Arguments

- `degrees::Bool`: Decide on which unit to return the angle in, degrees or radians. (default
                   is radians).

# Returns

- Pitch angle of the particle, either in radians (default) or degrees.
"""
function pitchangle(v, B; degrees::Bool=false)
    v = collect(v)
    B = collect(B)

    v_mag = norm(v)
    B_mag = norm(B)

    if v_mag ≤ eps(v_mag) || B_mag ≤ eps(B_mag)
        return 0.0
    end

    cross_vB = cross(v, B)
    cross_vB_mag = norm(cross_vB)
    dot_vB = dot(v, B)

    α = atan(cross_vB_mag, dot_vB)

    return degrees ? rad2deg(α) : α
end


"""
    get_v0(magnetic_field, r0, E_eV, m; α_frac=1)

Calculate the velocity vector for a particle in a magnetic field with some pitch angle.

A function that takes in a magnetic field model, a initial position, energy and mass, to
calculate an initial velocity that is within the loss cone, meaning that the particle will
precipitate.

# Arguments

- `magnetic_field`: Function describing the magnetic field, taking three positional values
                    as the argument.
- `r_0`: The initial position of the particle [m].
- `E_eV`: Energy of the particle [eV].
- `m`: Mass of the particle [kg].

# Keyword Arguments

- `α_frac`: Factor multiplied with the loss cone boundary.
"""
function get_v0(magnetic_field, r0, E_eV; α_frac=1)

    α_frac ≥ 0 || throw(ArgumentError("α_frac must be nonnegative"))

    if α_frac > 1
        @warn("The particle is now outside the loss cone.")
    end

    α = losscone_angle(magnetic_field, r0) * α_frac
    vy0 = sin(α) * velocity_from_kinetic_energy(E_eV, mₑ)
    vz0 = cos(α) * velocity_from_kinetic_energy(E_eV, mₑ)
    v0 = (0.0, vy0, vz0)

    return v0
end
