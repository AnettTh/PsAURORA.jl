"""
    parallel_velocity(particle::ParticleState, λ_grid::AbstractVector)

Calculate parallel velocity as a function of latitude in a magnetic field.

Uses `ParticleState` do get the parallel velocity for a given particle on a latitude-grid
`λ_grid`.

# Arguments

- `particle`: Structure holding the particle state.
- `λ_grid`: Magnetic latitude grid [rad].

# Returns

- A grid of parallel velocities corresponding to the latitude grid `λ_grid`.
"""
function parallel_velocity(particle::ParticleState, λ_grid::AbstractVector)

    v_parallel = zeros(length(λ_grid))

    # Find parallel velocity for each position
    for (i, λ) in enumerate(λ_grid)
        B_λ = particle.magnetic_field(particle.L, λ)
        B_λ_mag = norm(B_λ)

        α_λ = pitch_angle_at_λ(particle.α_eq, particle.B_eq, B_λ_mag)
        v_parallel[i] = particle.v * cos(α_λ)
    end

    return v_parallel
end


"""
    old_parallel_velocity(v, B)

Calculate the component of the velocity parallel to the magnetic field.

# Arguments

- `v`: Velocity vector.
- `B`: Magnetic-field vector.

# Returns

- The velocity vector parallel to `B`.

# Throws
- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function old_parallel_velocity(v, B)
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
function old_perpendicular_velocity(v, B)
    return v - old_parallel_velocity(v, B)
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
function old_parallel_speed(v, B)
    return norm(old_parallel_velocity(v, B))
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
function old_perpendicular_speed(v, B)
    return norm(old_perpendicular_velocity(v, B))
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
