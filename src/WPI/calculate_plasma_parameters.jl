using AURORA; dipole_field
using LinearAlgebra: norm, dot, cross
using StaticArrays

# TODO: update doc-string
"""
    losscone_angle(L; degree::Bool=false)
    losscone_angle(
    magnetic_field::AbstractMagneticField,
    L,
    ϕ0;
    degrees::Bool=false
    )

Calculate the equatorial pitch angle for the loss cone.

Using either the dipole approximation, or an arbitrary magnetic field funciton, as well as
the equatorial position, the function calculates the angle of the loss cone, i.e. the pitch
angle a particle in this position will need to have for it to precipitate into the
ionosphere.

# Arguments

- `L`: L-shell for where to find the loss cone.
- `magnetic_field`: Structure holding the magnetic field, doing the calculations without
  using the dipole approximation.
- `ϕ0`: Longitude of the field-line [rad].

# Keyword Arguments

- `degrees::Bool`: Choose units for the losscone angle, default is `false` (i.e. radians)
  and passing `true` will return output in degrees.

# Returns

- Loss cone angle for the given magnetic field and equatorial distance.

# Throws

- `ArgumentError`: If the equatorial distance is defined wrong.
- `ArgumentError`: If `B_ratio` is non-physical.
"""
function losscone_angle(
    fieldline::AbstractFieldLine;
    degrees::Bool=false,
    hemisphere::Symbol=:both
)
    # Equatorial fieldstrength already exist
    B_eq = fieldline.B_eq

    # Find the fieldstrength at the mirrorpoint
    B_mirror = hemisphere == :north ? fieldline.Bs[end] :
               hemisphere == :south ? fieldline.Bs[1] :
               hemisphere == :both  ? maximum(fieldline.Bs) :
               throw(ArgumentError("hemisphere must be :north, :south or :both"))

    B_ratio = B_eq / B_mirror

    if !(0.0 < B_ratio ≤ 1.0)
        throw(ArgumentError("Non-valid B_eq/B_mirror ratio: $B_ratio"))
    end

    # Calculate the losscone angle (TODO: Add reference for formula)
    α_rad = asin(sqrt(B_ratio))
    return degrees ? rad2deg(α_rad) : α_rad
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

- `ArgumentError`: If the given mass is zero or if the particles energy is negative.
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
    get_v0_from_Eμ(
    magnetic_field,
    r0,
    E_eV,
    μ0;
    ϕ=0.0,
    towards_equator::Bool=true,
    relativistic::Bool=false
)

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
- `μ0`: The cosine of the pitch-angle `α` of the particle, i.e. some value between -1 and 1.

# Keyword Arguments

- `ϕ`: The phase of the gyration, default is 0.0.
- `towards_equator`: Decides if the valocity is parallel or anti-parallel with the magnetic
  field, parallel being towards the equator IF in the northern hemisphere. The default is
  `true`.
- `relativistic`: Option to correct for relativistic effects, default is `false`.
"""
function get_v0_from_Eμ(
    magnetic_field,
    r0,
    E_eV,
    μ0;
    ϕ=0.0,
    towards_equator::Bool=true,
    relativistic::Bool=false
)

    -1 ≤ μ0 ≤ 1 || throw(ArgumentError("μ0 must be between -1 and 1"))

    v_mag = velocity_from_kinetic_energy(E_eV, mₑ; relativistic=relativistic)

    # Construct the magnetic basis
    B = magnetic_field(r0...)

    b, e1, e2, _ = magnetic_basis(B)

    # Adjust direction based on which hemisphere `r0` is in, meaning that if we are in the
    # northern hemisphere, towards equator means μ0 > 0
    if towards_equator
        μ0_corrected = r0[3] > 0 ? -abs(μ0) : abs(μ0)
    else
        μ0_corrected = r0[3] > 0 ? abs(μ0) : -abs(μ0)
    end

    # Determine parallel and perpendicular speeds
    v_parallel = μ0_corrected * v_mag
    v_perp = sqrt(1 - μ0_corrected^2) * v_mag

    # Initial velocity
    v0 = v_parallel .* b .+
         v_perp .* (cos(ϕ) .* e1 .+ sin(ϕ) .* e2)

    return Tuple(v0)
end
