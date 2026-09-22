using AURORA
using AURORA; z_ionosphere, mₑ, RE, BE

"""
    ParticleState{F<:Function}

Initial state of an electron for wave-particle interaction and time-of-flight calculations.

All derived quantities (speed, pitch angle, L-shell, loss cone angle, ionosphere latitude)
are computed from the initial position `r0`, pitch angle cosine `μ`, energy `E_eV` and
the magnetic field model.

# Fields

- `E_eV`: Kinetic energy [eV].
- `v`: Speed [m/s].
- `μ`: Pitch-angle cosine.
- `α0`: Pitch angle at initial position [rad].
- `r0`: Initial position in Cartesian coordinates [m].
- `B0`: Magnetic field vector at initial position [T].
- `b̂`: Unit vector along the magnetic field at initial position.
- `L`: L-shell number.
- `B_eq`: Magnetic field strength at the equator on this L-shell [T].
- `α_eq`: Pitch angle mapped to the equator via the adiabatic invariant [rad].
- `α_lc`: Loss cone angle at the initial position [rad].
- `λ_ionosphere`: Magnetic latitude where the L-shell intersects the ionosphere [rad].
- `magnetic_field`: Magnetic field model function `f(x, y, z)`.
"""
struct ParticleState{F<:Function}
    E_eV           :: Float64
    v              :: Float64
    μ              :: Float64
    α0             :: Float64
    r0             :: SVector{3, Float64}
    B0             :: SVector{3, Float64}
    b̂              :: SVector{3, Float64}
    L              :: Float64
    B_eq           :: Float64
    α_eq           :: Float64
    α_lc           :: Float64
    λ_ionosphere   :: Float64
    magnetic_field :: F
end


# TODO: Figure out how to better solve r0 and L, such that r0 can be anywhere and L is used for anything related to equator
"""
    ParticleState(E_eV, μ, r0, magnetic_field)

Construct a `ParticleState` from initial conditions.

Computes all derived quantities from the given energy, pitch angle cosine, initial position
and magnetic field model. The initial position `r0` is used to determine the L-shell and
all latitude-dependent quantities.

# Arguments

- `E_eV`: Energy of the electron [eV].
- `μ`: Cosine of the pitch angle, between -1 and 1.
- `r0`: Initial position vector in Cartesian coordinates [m].
- `magnetic_field`: Magnetic field function `f(x, y, z)`.

# Returns

- A fully initialized `ParticleState`.
"""
function ParticleState(E_eV, μ, r0, magnetic_field)

    v = velocity_from_kinetic_energy(E_eV, mₑ)
    B0 = magnetic_field(r0)
    b̂, _, _, _ = magnetic_basis(B0)
    r_mag = norm(r0)
    λ0 = asin(r0[3] / r_mag)
    L = r_mag / (RE * cos(λ0)^2)
    α_lc = losscone_angle(magnetic_field, r0)
    α0 = acos(abs(μ))
    B_eq = BE / L^3
    α_eq = pitch_angle_at_λ(α0, norm(B0), B_eq)
    λ_ionosphere = acos(sqrt((RE + z_ionosphere) / (L * RE))) # TODO: add throw

    return ParticleState(
        E_eV,
        v,
        μ,
        α0,
        SVector{3}(r0),
        SVector{3}(B0),
        SVector{3}(b̂),
        L,
        B_eq,
        α_eq,
        α_lc,
        λ_ionosphere,
        magnetic_field
    )
end
