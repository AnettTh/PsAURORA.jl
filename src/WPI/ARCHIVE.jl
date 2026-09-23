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

# Add this to ParticleState? Might be useful for boris-mover
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

# IDEA: add this to ParticleState?
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

- Gyrocenter position vector `SVector{3}` [m].

# Throws

- `ArgumentError`: Undefined if the magnitude of the magnetic field is zero.
"""
function gyrocenter(r, v, B, q, m)

    iszero(q) && throw(ArgumentError("Particle charge must be nonzero"))

    r_s = SVector{3}(r)
    v_s = SVector{3}(v)
    B_s = SVector{3}(B)

    B_mag = norm(B_s)
    iszero(B_mag) && throw(ArgumentError("Must have nonzero B"))

    b_hat = B_s / B_mag
    v_perp = v_s - dot(v_s, b_hat) * b_hat
    ρ = (m / (q * B_mag)) * cross(v_perp, b_hat)

    return r_s - ρ
end
