using AURORA
using LinearAlgebra: norm, dot, cross
using StaticArrays

# TODO: fix the 'old' functions to be multiple dispatch for struct, λ_grid
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
function gyro_frequency(B::AbstractVector, q, m)
    return gyro_frequency(norm(B), q, m)
end

function gyro_frequency(B_mag::Real, q, m)
    B_mag > 0 || throw(ArgumentError("Must have nonzero B"))
    return abs(q) * B_mag / m
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


# TODO: Add throws!
"""
    pitch_angle_at_λ(α_known, B_known, B_λ)

Calculate pitch angle for some magnetic field strength.

Given known pitch angle and magnetic field at some common location, the function calculates
the pitch angle corresponding to some magnetic field strength assuming the second adiabatic
invariant is conserved. If the particle mirrors, i.e. `sin²(α)` becomes larger than one,
`nothing` is returned.

# Arguments

- `α_known`: Known pitch angle at some location [rad].
- `B_known`: Known magnetic field strength at the same location as the pitch angle [T].
- `B_λ`: Magnetic field strength at latitude `λ` [T].

# Returns

- Pitch angle at latitude `λ` [rad].
"""
function pitch_angle_at_λ(α_known, B_known, B_λ)

    iszero(B_known) && throw(ArgumentError("B_known must be nonzero"))
    abs(α_known) ≥ 2π && throw(ArgumentError("α_known must be between ±2π"))

    sin2_α = sin(α_known)^2 * B_λ / B_known

    # Check for mirroring
    sin2_α > 1 && return nothing

    α = asin(sqrt(sin2_α))

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
# TODO: Add relativistic option
function velocity_from_kinetic_energy(E_eV, m)

    iszero(m) && throw(ArgumentError("Mass must be nonzero"))

    E_J = E_eV * eV_in_J

    ratio = E_J / (m * c₀^2)

    #if ratio ≥ 1
    #    throw(
    #        ArgumentError(
    #            "Kinetic energy is too large for non-relaticistic velocity approximation"
    #        )
    #    )
    #elseif ratio ≥ 0.1
    #    @warn "Relativistic corrections might be significant"
    #end

    m > 0 || throw(ArgumentError("Mass must be positive"))
    E_eV ≥ 0 || throw(ArgumentError("Kinetic energy must be nonnegative"))

    v = sqrt(2E_J / m)

    return v
end


"""
    get_v0_from_Eμ(magnetic_field, r0, E_eV, μ; ϕ=0.0, towards_equator::Bool=true)

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
- `μ`: The cosine of the pitch-angle `α` of the particle, i.e. some value between -1 and 1.

# Keyword Arguments

- `ϕ`: The phase of the gyration, default is 0.0.
- `towards_equator`: Decides if the valocity is parallel or anti-parallel with the magnetic
  field, parallel being towards the equator IF in the northern hemisphere. The default is
  `true`.
"""
function get_v0_from_Eμ(magnetic_field, r0, E_eV, μ; ϕ=0.0, towards_equator::Bool=true)

    -1 ≤ μ ≤ 1 || throw(ArgumentError("μ must be between -1 and 1"))

    v_mag = velocity_from_kinetic_energy(E_eV, mₑ)

    # Construct the magnetic basis
    B = magnetic_field(r0...)

    b, e1, e2, _ = magnetic_basis(B)

    # Adjust direction based on which hemisphere `r0` is in, meaning that if we are in the
    # northern hemisphere, towards equator means μ > 0
    if towards_equator
        μ_corrected = r0[3] > 0 ? -abs(μ) : abs(μ)
    else
        μ_corrected = r0[3] > 0 ? abs(μ) : -abs(μ)
    end

    # Determine parallel and perpendicular speeds
    v_parallel = μ_corrected * v_mag
    v_perp = sqrt(1 - μ_corrected^2) * v_mag

    # Initial velocity
    v0 = v_parallel .* b .+
         v_perp .* (cos(ϕ) .* e1 .+ sin(ϕ) .* e2)

    return Tuple(v0)
end

# TODO: Add and test whistler wave number with option for relativistic energies
# TODO: Add and test whistler wave dispersion relation with option for relativistic energies
# TODO: Add and test resonance condition with option for relativistic energies
