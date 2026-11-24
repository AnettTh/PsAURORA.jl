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
- `μ0`: Pitch-angle cosine.
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
- `ϕ0`: Longitude at equator [rad].
- `φ0`: Initial gyrophase of the particle [rad].
"""
struct ParticleState{F}
    E0_eV  :: Float64
    μ0     :: Float64
    α0     :: Float64
    r0     :: SVector{3, Float64}
    ϕ0     :: Float64
    φ0     :: Float64

    v0_mag :: Float64
    v0     :: SVector{3, Float64}
    B0     :: SVector{3, Float64}
    b̂      :: SVector{3, Float64}
    E0     :: SVector{3, Float64}

    r_eq    :: Float64
    B_eq    :: Float64
    α_eq    :: Float64
    α_eq_lc :: Float64

    magnetic_field :: F

    λ_ionosphere :: Float64

    γ :: Float64
end

"""
    ParticleState(
    E0_eV,
    μ0,
    fieldline::AbstractFieldLine;
    r0::Cartesian=fieldline.position_eq,
    relativistic::Bool=true,
    φ0=deg2rad(0.0),
)

Construct a `ParticleState` from initial conditions.

Computes all derived quantities from the given energy, pitch angle cosine, initial position
and magnetic field model, either with or without relativistic corrections. The initial
position `r0` is used to determine the L-shell and all latitude-dependent quantities.

# Arguments

- `E_eV`: Energy of the electron [eV].
- `μ0`: Cosine of the pitch angle, between -1 and 1.
- `r0`: Initial position vector in Cartesian coordinates [m].
- `magnetic_field`: Magnetic field function `f(x, y, z)`.
- `ϕ0`: The longitudinal position of the particle at the equator [rad].

# Keyword Arguments

- `relativistic`: Option to correct for relativistic effects, default is `false`.
- `φ0`: Initial gyrophase of the particle [rad].

# Returns

- A fully initialized `ParticleState`.
"""
function ParticleState(
    E0_eV,
    μ0,
    fieldline::AbstractFieldLine;
    r0::Cartesian=fieldline.position_eq,
    relativistic::Bool=true,
    φ0=deg2rad(0.0),
)
    magnetic_field = fieldline.magnetic_field

    # Check validity of initial position
    r_mag = norm((r0.x, r0.y, r0.z))
    r_mag ≤ RE && throw(ArgumentError("r0 must be outside the Earth's surface."))

    # Pitch angle and longitude at initial position
    α0 = acos(abs(μ0))
    ϕ0 = mod2pi(atan(r0.y, r0.x))

    # Get velocity and relativistic constant
    v0_mag, γ = velocity_from_kinetic_energy(
        E0_eV,
        mₑ;
        relativistic=relativistic,
        return_gamma=true,
    )

    # Fields at initial position
    B0 = magnetic_field(r0)
    b̂, e1, e2, _ = magnetic_basis(B0)

    # Add also an electric field, for consistency
    E0 = 0.0, 0.0, 0.0

    # Find initial velocity
    v_par = -v0_mag * μ0
    v_perp = v0_mag * sqrt(1 - μ0^2)
    v0 = SVector{3}(
        v_par .* b̂ .+
        v_perp .* cos(φ0) .* e1 .+
        v_perp .* sin(φ0) .* e2
    )

    # Values at equator, minimum(|B|), needed to fill the struct
    B_eq = fieldline.B_eq
    α_eq = pitch_angle(α0, norm(B0), B_eq) # Figure out if α_eq should be equal α0 here

    # Loss cone and ionospheric latitude for correct hemisphere
    hemisphere = μ0 < 0 ? :north : :south
    α_eq_lc = losscone_angle(fieldline; hemisphere=hemisphere)
    λ_ionosphere = abs(μ0 < 0 ? fieldline.λs[end] : fieldline.λs[1])

    return ParticleState(
        E0_eV,
        μ0,
        α0,
        SVector(r0.x, r0.y, r0.z),
        ϕ0,
        φ0,

        v0_mag,
        v0,
        SVector{3}(B0),
        SVector{3}(b̂),
        SVector{3}(E0),

        fieldline.r_eq,
        B_eq,
        α_eq,
        α_eq_lc,

        magnetic_field,

        λ_ionosphere,

        γ,
    )
end


function Base.show(io::IO, p::ParticleState)
    println(io, "ParticleState:")
    println(io, "  E       = $(round(p.E0_eV/1e3, digits=2)) keV")
    println(io, "  α₀      = $(round(rad2deg(p.α0), digits=2))°")
    println(io, "  α_eq    = $(round(rad2deg(p.α_eq), digits=5))°")
    println(io, "  α_lc    = $(round(rad2deg(p.α_eq_lc), digits=5))°")
    println(io, "  R_eq    = $(round(p.r_eq/RE, digits=2)) RE")
    println(io, "  ϕ0      = $(round(rad2deg(p.ϕ0), digits=1))° " *
                "($(round(longitude_to_MLT(p.ϕ0, degrees=false), digits=1)) MLT)")
    println(io, "  γ       = $(round(p.γ, digits=4))")
    println(io, "  B_field = $(typeof(p.magnetic_field))")
end
