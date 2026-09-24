using AURORA

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
