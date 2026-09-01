using AURORA
using LinearAlgebra: norm

"""
    gyro_frequency(B, q, m)

Calculate the gyro-frequency of a test particle in a magnetic field.

# Arguments

- `B`: Magnetic field vector `[Bx, By, Bz]` [T].
- `q`: The charge of the test particle [C].
- `m`: The mass of the test particle [kg].

# Returns

- Gyrofrequency of the test particle [rad/s]

# Throws

- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function gyro_frequency(B, q, m)
    B_mag = norm(B)

    iszero(B_mag) && throw(ArgumentError("Must have nonzero B"))

    return abs(q) * B_mag / m
end


"""
    parallel_velocity(v, B)

Calculate the component of the velocity parallel to the magnetic field.

# Arguments

- `v`: Velocity vector.
- `B`: Magnetic-field vector.

# Returns

- The velocity vector parallel to `B`.

# Throws
- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function parallel_velocity(v, B)
    B_mag2 = dot(B, B)

    iszero(B_mag2) && throw(ArgumentError("Must have nonzero B"))

    return (dot(v, B) / B_mag2) * B
end


"""
    perpendicular_velocity(v, B)

Calculate the component of the velocity perpendicular to the magnetic field.

# Arguments

- `v`: Velocity vector.
- `B`: Magnetic-field vector.

# Returns

- The velocity vector perpendicular to `B`.

# Throws
- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function perpendicular_velocity(v, B)
    return v - parallel_velocity(v, B)
end


"""
    parallel_speed(v, B)

Calculate the magnitude of the velocity component parallel to the magnetic field.

# Arguments

- `v`: Velocity vector.
- `B`: Magnetic-field vector.

# Returns

- The magnitude of the velocity parallel to `B`.

# Throws
- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function parallel_speed(v, B)
    return norm(parallel_velocity(v, B))
end


"""
    perpendicular_speed(v, B)

Calculate the magnitude of the velocity component perpendicular to the magnetic field.

# Arguments

- `v`: Velocity vector.
- `B`: Magnetic-field vector.

# Returns

- The magnitude of the velocity perpendicular to `B`.

# Throws
- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function perpendicular_speed(v, B)
    return norm(perpendicular_velocity(v, B))
end


"""
    larmor_radius(m, q, v, B)

Calculate the larmor radius of a test particle in a magnetic field.

# Arguments

- `m`: The mass of the test particle [kg].
- `q`: The charge of the test particle [C].
- `v`: Velocity vector of the test particle `[vx, vy, vz]` [m/s].
- `B`: Magnetic field vector `[Bx, By, Bz]` [T].

# Returns

- Larmor radius of the test particle [m]

# Throws

- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function larmor_radius(m, q, v, B)
    B_mag = norm(B)
    B_mag ≤ eps(B_mag) && throw(ArgumentError("Must have nonzero B"))
    iszero(q) && throw(ArgumentError("Particle charge must be nonzero"))

    v_perp = perpendicular_speed(v, B)

    return (m * v_perp) / (abs(q) * B_mag)
end


# IDEA: Can be made more efficient by staticarrays.jl
"""
    gyrocenter(r, v, B, q, m)

Calculate the position of the gyrocenter of a test particle in a magnetic field.

# Arguments

- `r`: Position vector of the test particle `[x, y, z]` [m].
- `v`: Velocity vector of the test particle `[vx, vy, vz]` [m/s].
- `B`: Magnetic field vector `[Bx, By, Bz]` [T].
- `q`: The charge of the test particle [C].
- `m`: The mass of the test particle [kg].

# Returns

- Gyrocenter position vector `[x, y, z]` [m].

# Throws

- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function gyrocenter(r, v, B, q, m)

    # If tuples, convert to arrays to allow for linear algebra
    r = collect(r)
    v = collect(v)
    B = collect(B)

    B_mag = norm(B)

    iszero(B_mag) && throw(ArgumentError("Must have nonzero B"))
    iszero(q) && throw(ArgumentError("Particle charge must be nonzero"))

    b_hat = B / B_mag
    v_perp = v - dot(v, b_hat) * b_hat
    ρ = (m / (q * B_mag)) * cross(v_perp, b_hat)

    return r - ρ
end

# TODO: Change input to L-shell?
"""
    losscone_angle(dipole_field, r_eq; degree::Bool=false)

Calculate the equatorial pitch angle for the loss cone.

Using the defined magnetic field and the equatorial position, the function calculates the
angle of the loss cone, i.e. the pitch angle a particle in this position will need to have
for it to precipitate into the ionosphere.

# Arguments

- `dipole_field`: Function describing the magnetic field, taking three positional values
                    as the argument.
- `r_eq`: The position in the equatorial plane for where to find the loss cone [m].

# Keyword Arguments

- `degrees::Bool`: Choose units for the pitch angle, default is `false` (i.e. radians) and
                   passing `true` will allow for input in degrees.

# Returns

- Loss cone angle for the given magnetic field and equatorial distance.

# Throws

- `ArgumentError`: If the equatorial distance is defined wrong.
"""
function losscone_angle(dipole_field, r_eq; degrees::Bool=false)
    r_mirror = RE + z_ionosphere
    r_eq_mag = norm(r_eq)

    iszero(r_eq_mag) && throw(ArgumentError("Equatorial distance `r_eq` must be nonzero"))

    mirror_arg = r_mirror / r_eq_mag

    λ_mirror = acos(sqrt(mirror_arg))

    x_mirror = r_mirror * cos(λ_mirror)
    y_mirror = 0.0
    z_mirror = r_mirror * sin(λ_mirror)

    B_eq = norm(dipole_field(r_eq...))
    B_mirror = norm(dipole_field(x_mirror, y_mirror, z_mirror))

    B_ratio = B_eq / B_mirror

    if !(0.0 < B_ratio ≤ 1.0)
        throw(ArgumentError("The chosen positions return non-valid magnetic field vectors"))
    end

    α_rad = asin(sqrt(B_ratio))

    α = degrees ? rad2deg(α_rad) : α_rad

    return α
end



"""
    velocity_from_kinetic_energy(E_eV, m)

Calculate magnitude of velocity vector from kinetic energy.

Takes in kinetic energy given as electron volts and mass of some particle, converts the
energy to Joules and then finds the magnitude of the velocity of the given particle.

# Arguments

- `E_eV`: Kinetic energy of the particle [eV].
- `m`: Mass of the particle [kg].

# Returns

- Magnitude of the velocity of the particle [m/s].

# Throws

- `ArgumentError`: If the given mass is zero or if the particles speed is faster than the
                   speed of light.
"""
function velocity_from_kinetic_energy(E_eV, m)

    iszero(m) && throw(ArgumentError("Mass must be nonzero"))

    E_J = E_eV * eV_in_J

    ratio = E_J / (m * c^2)

    if ratio ≥ 1
        throw(
            ArgumentError(
                "Kinetic energy is too large for non-relaticistic velocity approximation"
            )
        )
    elseif ratio ≥ 0.1
        @warn "Relativistic corrections might be significant"
    end

    m > 0 || throw(ArgumentError("Mass must be positive"))
    E_eV ≥ 0 || throw(ArgumentError("Kinetic energy must be nonnegative"))

    v = sqrt(2E_J / m)

    return v
end


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
    quarter_bounceperiod(L, E_eV, m, θ; degrees::Bool=false)

Calculate the quarter bounce period for a particle in a dipole magnetic field.

Using an approximation for the change in arclength with respect to change in latitude for a
mirroring particle, the quarter of one full bounce period can be calculated. The function
takes in the equatorial position, i.e. the L-shell, as well as the energy of the particle,
the mass of the particle and the pitch angle, and returns an approximation for how long it
takes for the particle to travel along a dipolar magnetic field and to the mirror-point.

# Arguments

- `L`: Initial distance of the particle, given as the L-shell number.
- `E_eV`: Energy of the particle [eV].
- `m`: Mass of the particle [kg].
- `θ`: Pitch-angle of the particle, either in radians or degrees.

# Keyword Arguments

- `degrees::Bool`: Choose units for the pitch angle, default is `false` (i.e. radians) and
                   passing `true` will allow for input in degrees.
"""
function quarter_bounceperiod(L, E_eV, m, θ; degrees::Bool=false)

    θ_rad = degrees ? deg2rad(θ) : θ

    v = velocity_from_kinetic_energy(E_eV, m)

    # Gamma-function used as an approximation to a unsolvable integral
    Γ = 1.30 - 0.56 * sin(θ_rad)

    return ((L * RE) / v) * Γ
end


"""
    average_driftvelocity(L, E_eV, q, θ)

Calculate the average drift velocity for a particle in a dipole magnetic field.

Using an approximation for the change in angular velocity, integrated over a quarter of
a bounce period allows us to find the average drift velocity. The function takes in the
equatorial position, i.e. the L-shell, as well as the energy of the particle, the charge of
the particle and the pitch angle, and returns an approximation for the drift velocity of the
particle.

# Arguments

- `L`: Initial distance of the particle, given as the L-shell number.
- `E_eV`: Energy of the particle [eV].
- `q`: Charge of the particle [C].
- `θ`: Pitch-angle of the particle, either in radians or degrees.

# Keyword Arguments

- `degrees::Bool`: Choose units for the pitch angle, default is `false` (i.e. radians) and
                   passing `true` will allow for input in degrees.
"""
function average_driftvelocity(L, E_eV, q, θ; degrees::Bool=false)

    θ_rad = degrees ? deg2rad(θ) : θ

    E_J = E_eV * eV_in_J

    # Gamma-function used as an approximation to a unsolvable integral
    Γ = 0.35 - 0.15 * sin(θ_rad)

    return (3 * L^2 * E_J * Γ) / (2 * q * BE * RE)
end


"""
    total_drift(L, E_eV, q, m, θ; degrees::Bool=false)

Calculate total drift based on bounce period and average drift velocity.

Using an approximation for change in arclength with respect to change in latitude and the
change in angular velocity, integrated over a quarter of a bounce period allows us to find
the total drift of a particla. The function takes in the equatorial position, i.e. the
L-shell, as well as the energy of the particle, the mass and charge of the particle and the
pitch angle, and returns an approximation for the total drift of the particle.

# Arguments

- `L`: Initial distance of the particle, given as the L-shell number.
- `E_eV`: Energy of the particle [eV].
- `q`: Charge of the particle [C].
- `m`: Mass of the particle [kg].
- `θ`: Pitch-angle of the particle, either in radians or degrees.

# Keyword Arguments

- `degrees::Bool`: Choose units for the pitch angle, default is `false` (i.e. radians) and
                   passing `true` will allow for input in degrees.

"""
function total_drift(L, E_eV, q, m, θ; degrees::Bool=false)

    τ = quarter_bounceperiod(L, E_eV, m, θ; degrees=degrees)
    v_d = average_driftvelocity(L, E_eV, q, θ; degrees=degrees)

    return τ * v_d
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
function get_v0(magnetic_field, r0, E_eV, m; α_frac=1)

    α_frac ≥ 0 || throw(ArgumentError("α_frac must be nonnegative"))

    if α_frac > 1
        @warn("The particle is now outside the loss cone.")
    end

    α = losscone_angle(magnetic_field, r0) * α_frac
    vy0 = sin(α) * velocity_from_kinetic_energy(E_eV, m)
    vz0 = cos(α) * velocity_from_kinetic_energy(E_eV, m)
    v0 = (0.0, vy0, vz0)

    return v0
end
