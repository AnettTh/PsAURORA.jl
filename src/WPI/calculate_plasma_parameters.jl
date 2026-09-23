using AURORA; dipole_field
using LinearAlgebra: norm, dot, cross
using StaticArrays

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
# NOTE: This can be archived if the boris-mover in a simple way can use L, λ for B
function gyro_frequency(B::AbstractVector, q, m)
    return gyro_frequency(norm(B), q, m)
end

function gyro_frequency(B_mag::Real, q, m)
    B_mag > 0 || throw(ArgumentError("Must have nonzero B"))
    return abs(q) * B_mag / m
end


"""
    losscone_angle(L; degree::Bool=false)

Calculate the equatorial pitch angle for the loss cone.

Using the defined magnetic field and the equatorial position, the function calculates the
angle of the loss cone, i.e. the pitch angle a particle in this position will need to have
for it to precipitate into the ionosphere.

# Arguments

- `dipole_field`: Function describing the magnetic field, taking three positional values
  as the argument.
- `L`: L-shell for where to find the loss cone.

# Keyword Arguments

- `degrees::Bool`: Choose units for the losscone angle, default is `false` (i.e. radians)
  and passing `true` will return output in degrees.

# Returns

- Loss cone angle for the given magnetic field and equatorial distance.

# Throws

- `ArgumentError`: If the equatorial distance is defined wrong.
"""
function losscone_angle(L; degrees::Bool=false)
    r_mirror = RE + z_ionosphere
    r_eq_mag = L * RE

    iszero(r_eq_mag) && throw(ArgumentError("Equatorial distance `r_eq` must be nonzero"))

    mirror_arg = r_mirror / r_eq_mag

    λ_mirror = acos(sqrt(mirror_arg))

    x_mirror = r_mirror * cos(λ_mirror)
    y_mirror = 0.0
    z_mirror = r_mirror * sin(λ_mirror)

    B_eq = norm(dipole_field(r_eq_mag, 0.0, 0.0))
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
    velocity_from_kinetic_energy(
    E_eV,
    m;
    relativistic::Bool=false,
    return_gamma::Bool=false
    )

Calculate magnitude of velocity vector from kinetic energy.

Takes in kinetic energy given as electron volts and mass of some particle, converts the
energy to Joules and then finds the magnitude of the velocity of the given particle, either
with or without relativistic corrections.

# Arguments

- `E_eV`: Kinetic energy of the particle [eV].
- `m`: Mass of the particle [kg].

# Keyword Arguments

- `relativistic`: Option to correct for relativistic effects, default is `false`.
- `return_gamma`: Option to also return the relativistic factor `γ`, default is `false`.

# Returns

- Magnitude of the velocity of the particle [m/s].

# Throws

- `ArgumentError`: If the given mass is zero or if the particles speed is faster than the
  speed of light.
"""
function velocity_from_kinetic_energy(
    E_eV,
    m;
    relativistic::Bool=false,
    return_gamma::Bool=false
    )

    m > 0 || throw(ArgumentError("Mass must be positive"))
    E_eV ≥ 0 || throw(ArgumentError("Kinetic energy must be nonnegative"))

    E_J = E_eV * eV_in_J

    if relativistic
        γ = E_J / (m * c₀^2) + 1
        v = c₀ * sqrt(1 - 1/γ^2)
    else
        # NOTE: Add this to the code later? For testing, it's fine to leave it out
        #ratio = E_J / (m * c₀^2)
        #ratio ≥ 1 && throw(ArgumentError(
        #    "Kinetic energy too large for non-relativistic approximation"
        #))
        #ratio ≥ 0.1 && @warn "Relativistic corrections might be significant"

        v = sqrt(2E_J / m)
        γ = 1.0
    end

    if return_gamma
        return v, γ
    else
        return v
    end
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
- `relativistic`: Option to correct for relativistic effects, default is false.
"""
function get_v0_from_Eμ(magnetic_field, r0, E_eV, μ; ϕ=0.0, towards_equator::Bool=true, relativistic::Bool=false)

    -1 ≤ μ ≤ 1 || throw(ArgumentError("μ must be between -1 and 1"))

    v_mag = velocity_from_kinetic_energy(E_eV, mₑ; relativistic=relativistic)

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
