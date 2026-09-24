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

#===================================Run preferred method===================================#
# IDEA: Make this into 'AbstractPropagation' with WPI also!!
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
