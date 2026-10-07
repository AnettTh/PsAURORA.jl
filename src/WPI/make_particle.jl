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
- `γ`: Lorentz factor.
- `ϕ`: Longitude [rad].
"""
struct ParticleState{F}
    E_eV           :: Float64
    v              :: Float64
    μ              :: Float64
    α0             :: Float64
    r0             :: SVector{3, Float64}   # TODO: Change this to Cartesian, and maybe add Spherical?
    v0             :: SVector{3, Float64}
    B0             :: SVector{3, Float64}
    b̂              :: SVector{3, Float64}
    L              :: Float64
    B_eq           :: Float64
    α_eq           :: Float64
    α_lc           :: Float64
    λ_ionosphere   :: Float64
    magnetic_field :: F             # TODO: See if this can be done without storing the magnetic field here also
    B_interpolated :: Any
    s0             :: Float64
    γ              :: Float64
    ϕ              :: Float64
    ϕ_eq           :: Float64
    φ0             :: Float64
end

function find_R_max(particle::ParticleState; kwargs...)
    return find_R_max(
        particle.magnetic_field,
        particle.L,
        0.0,
        particle.ϕ_eq;
        kwargs...
    )
end

# TODO: Test this for non-equatorial r0, check that values makes sense physically
# TODO: See if r0 actually is particle r0, or if it just represents the equatorial distance of the field-line in question
# TODO: Separate between L and Spherical, need r instead!
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
- `ϕ_eq`

# Keyword Arguments

- `precomputed_magnetic_field`: #TODO: describe
- `relativistic`: Option to correct for relativistic effects, default is false.
- `φ`: #TODO: describe

# Returns

- A fully initialized `ParticleState`.
"""
function ParticleState(
    E_eV,
    μ,
    r0::Cartesian,
    magnetic_field::AbstractMagneticField,
    ϕ_eq;
    precomputed_magnetic_field::Bool=false,
    relativistic::Bool=false,
    φ0=deg2rad(0.0)
)

    # Check validity of initial position
    r_mag = norm((r0.x, r0.y, r0.z))
    r_mag ≤ RE && throw(ArgumentError("r0 must be outside the Earth's surface."))

    # Get velocity and relativistic constant
    v, γ = velocity_from_kinetic_energy(
        E_eV,
        mₑ;
        relativistic=relativistic,
        return_gamma=true
    )

    # Specify magnetic field
    B0 = magnetic_field(r0)
    b̂, e1, e2, _ = magnetic_basis(B0)

    # Find initial velocity
    v_par = v * μ
    v_perp = v * sqrt(max(1 - μ^2, 0.0))
    v0 = SVector{3}(        # TODO: This has some issue, figure out what sqrt of negative number
        v_par .* b̂ .+
        v_perp .* cos(φ0) .* e1 .+
        v_perp .* sin(φ0) .* e2
    )

    # Latitude of initial position and L-shell
    λ0 = asin(r0.z / r_mag)
    L = r_mag / (RE * cos(λ0)^2)

    # TODO: test this!!
    ϕ = atan(r0.z, r0.x)

    # Loss-cone angle at the equator
    # IDEA: Look into the definition of α from Hsieh 2022
    if magnetic_field isa DipoleMagneticField
        α_lc = losscone_angle(L)        # TODO: Test this for off-equatorial positions, to check for bugs! It should be larger than α_eq?
    else
        α_lc = NaN
        # IDEA: Include also non-dipolar losscone, if not too expensive? Check what it is used for first
    end

    # Current position pitch-angle and equatorial pitch-angle
    α0 = acos(abs(μ))
    B_eq = BE / L^3
    α_eq = pitch_angle_at_λ(α0, norm(B0), B_eq)

    # Latitude of the ionosphere based on the L-shell
    λ_ionosphere = acos(sqrt((RE + z_ionosphere) / (L * RE)))

    # Pre-compute magnetic field
    if precomputed_magnetic_field
        B_interpolated, s_grid = get_magnetic_field(magnetic_field, L, ϕ_eq)
        s0 = s_grid[argmin([norm(B_interpolated(s)) for s in s_grid])]
    else
        B_interpolated = nothing
        s0 = 0.0
    end

    return ParticleState(
        E_eV,
        v,
        μ,
        α0,
        SVector(r0.x, r0.y, r0.z),
        v0,
        SVector{3}(B0),
        SVector{3}(b̂),
        L,
        B_eq,
        α_eq,
        α_lc,
        λ_ionosphere,
        magnetic_field,
        B_interpolated,
        s0,
        γ,
        ϕ,
        ϕ_eq,
        φ0
    )
end

# TODO: Fix this to something useful
function Base.show(io::IO, p::ParticleState)
    println(io, "ParticleState:")
    println(io, "  E        = $(round(p.E_eV/1e3, digits=2)) keV")
    println(io, "  α₀       = $(round(rad2deg(p.α0), digits=2))°")
    println(io, "  α_eq     = $(round(rad2deg(p.α_eq), digits=2))°")
    println(io, "  α_lc     = $(round(rad2deg(p.α_lc), digits=2))°")
    println(io, "  L        = $(round(p.L, digits=2))")
    println(io, "  ϕ_eq     = $(round(rad2deg(p.ϕ_eq), digits=1))° ($(round(longitude_to_MLT(p.ϕ_eq, degrees=false), digits=1)) MLT)")
    println(io, "  γ        = $(round(p.γ, digits=4))")
    println(io, "  B_field  = $(typeof(p.magnetic_field))")
    print(io,   "  B_interp = $(isnothing(p.B_interpolated) ? "not precomputed" : "precomputed")")
end
