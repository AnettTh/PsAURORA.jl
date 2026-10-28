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
- `v0`: Initial velocity vector [m/s].
- `B0`: Magnetic field vector at initial position [T].
- `E0`: Electric field vector at initial position [V/m].
- `b̂`: Unit vector along the magnetic field at initial position.
- `L`: L-shell number.
- `B_eq`: Magnetic field strength at the equator on this L-shell [T].
- `α_eq`: Pitch angle mapped to the equator via the adiabatic invariant [rad].
- `α_lc`: Loss cone angle at the initial position [rad].
- `λ_ionosphere`: Magnetic latitude where the L-shell intersects the ionosphere [rad].
- `magnetic_field`: Magnetic field model function `f(x, y, z)`.
- `γ`: Lorentz factor.
- `ϕ_eq`: Longitude at equator [rad].
- `φ0`: Initial gyrophase of the particle [rad].
"""
struct ParticleState{F}
    E_eV           :: Float64
    v              :: Float64
    μ              :: Float64
    α0             :: Float64
    r0             :: SVector{3, Float64}
    v0             :: SVector{3, Float64}
    B0             :: SVector{3, Float64}
    E0             :: SVector{3, Float64}
    b̂              :: SVector{3, Float64}
    L              :: Float64
    B_eq           :: Float64
    α_eq           :: Float64
    α_lc           :: Float64
    λ_ionosphere   :: Float64
    magnetic_field :: F
    γ              :: Float64
    ϕ_eq           :: Float64
    φ0             :: Float64
end

# TODO: Check if this is necessary anymore
function find_R_max(particle::ParticleState; kwargs...)
    return find_R_max(
        particle.magnetic_field,
        Spherical(particle.L, 0.0, particle.ϕ_eq);
        kwargs...
    )
end

# TODO: Test this for non-equatorial r0, check that values makes sense physically
# TODO: See if r0 actually is particle r0, or if it just represents the equatorial distance of the field-line in question
# TODO: Look into if the particle is defined correctly also for a tsyganenko field
"""
    ParticleState(E_eV, μ, r0, magnetic_field, ϕ_eq; relativistic::Bool=false)

Construct a `ParticleState` from initial conditions.

Computes all derived quantities from the given energy, pitch angle cosine, initial position
and magnetic field model, either with or without relativistic corrections. The initial
position `r0` is used to determine the L-shell and all latitude-dependent quantities.

# Arguments

- `E_eV`: Energy of the electron [eV].
- `μ`: Cosine of the pitch angle, between -1 and 1.
- `r0`: Initial position vector in Cartesian coordinates [m].
- `magnetic_field`: Magnetic field function `f(x, y, z)`.
- `ϕ_eq`: The longitudinal position of the particle at the equator [rad].

# Keyword Arguments

- `relativistic`: Option to correct for relativistic effects, default is `false`.
- `φ0`: Initial gyrophase of the particle [rad].

# Returns

- A fully initialized `ParticleState`.
"""
function ParticleState(
    E_eV,
    μ,
    r0::Cartesian,
    magnetic_field::AbstractMagneticField,
    ϕ_eq;
    relativistic::Bool=false,
    φ0=deg2rad(0.0),
)

    # Check validity of initial position
    r_mag = norm((r0.x, r0.y, r0.z))
    r_mag ≤ RE && throw(ArgumentError("r0 must be outside the Earth's surface."))

    # Get velocity and relativistic constant
    v, γ = velocity_from_kinetic_energy(
        E_eV,
        mₑ;
        relativistic=relativistic,
        return_gamma=true,
    )

    # Specify magnetic field
    B0 = magnetic_field(r0)
    b̂, e1, e2, _ = magnetic_basis(B0)

    # Add also an electric field, for consistency
    E0 = 0.0, 0.0, 0.0

    # Find initial velocity
    v_par = -v * μ
    v_perp = v * sqrt(1 - μ^2)
    v0 = SVector{3}(
        v_par .* b̂ .+
        v_perp .* cos(φ0) .* e1 .+
        v_perp .* sin(φ0) .* e2
    )

    # Find point of max R to use as L-shell
    L = find_R_max(magnetic_field, r0) / RE

    # Loss-cone angle at the equator of this L-shell
    α_lc = losscone_angle(magnetic_field, L, ϕ_eq)

    # Current position pitch-angle and equatorial pitch-angle
    α0 = acos(abs(μ))

    # Latitude of the ionosphere based on the L-shell
    λ_ionosphere = acos(sqrt((RE + z_ionosphere) / (L * RE)))

    # Find field-strength and pitch angle at the equator
    # TODO: These are not correct for the tsyganenko-model, look into how to handle it
    B_eq = BE / L^3
    α_eq = pitch_angle_at_λ(α0, norm(B0), B_eq)

    return ParticleState(
        E_eV,
        v,
        μ,
        α0,
        SVector(r0.x, r0.y, r0.z),
        v0,
        SVector{3}(B0),
        SVector{3}(E0),
        SVector{3}(b̂),
        L,
        B_eq,
        α_eq,
        α_lc,
        λ_ionosphere,
        magnetic_field,
        γ,
        ϕ_eq,
        φ0
    )
end

function ParticleState(
    E_eV,
    μ,
    r0::Spherical,
    magnetic_field::AbstractMagneticField,
    ϕ_eq;
    kwargs...
)
    return ParticleState(E_eV, μ, Cartesian(r0), magnetic_field, ϕ_eq; kwargs...)
end

function Base.show(io::IO, p::ParticleState)
    println(io, "ParticleState:")
    println(io, "  E        = $(round(p.E_eV/1e3, digits=2)) keV")
    println(io, "  α₀       = $(round(rad2deg(p.α0), digits=2))°")
    println(io, "  α_eq     = $(round(rad2deg(p.α_eq), digits=2))°")
    println(io, "  α_lc     = $(round(rad2deg(p.α_lc), digits=2))°")
    println(io, "  L        = $(round(p.L, digits=2))")
    println(io, "  ϕ_eq     = $(round(rad2deg(p.ϕ_eq), digits=1))° " *
                "($(round(longitude_to_MLT(p.ϕ_eq, degrees=false), digits=1)) MLT)")
    println(io, "  γ        = $(round(p.γ, digits=4))")
    println(io, "  B_field  = $(typeof(p.magnetic_field))")
end
