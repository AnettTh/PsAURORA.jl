using Dates
using TsyganenkoModels
using StaticArrays
using LinearAlgebra

"""
    magnetic_basis(B)

Construct an orthonormal magnetic-field basis from a magnetic-field vector.

The first basis vector, `̂b`, is parallel to the magnetic field, while `e1` and `e2` spans
the plane perpendicular to the magnetic field.
# Arguments

- `B`: Magnetic-field vector in Cartesian coordinates.

# Returns

A tuple `(b̂, e1, e2, B_mag)` containing:

- `b̂`: Unit vector parallel to the magnetic field.
- `e1`: Unit vector perpendicular to `b̂`.
- `e2`: Unit vector perpendicular to both `b̂` and `e1`.
- `B_mag`: Magnitude of the magnetic-field vector.

# Throws

- `ArgumentError`: If the magnetic-field magnitude is zero.
"""
function magnetic_basis(B)
    B = collect(B)

    B_mag = norm(B)
    B_mag > eps() || throw(ArgumentError("Magnetic-field magnitude must be non-zero"))

    b̂ = B ./ B_mag

    # Use x-direction as one vector if b̂ is not too close, else use y-direction
    if abs(b̂[1]) < 0.9
        safe_vector = [1.0, 0.0, 0.0]
    else
        safe_vector = [0.0, 1.0, 0.0]
    end


    e1 = cross(b̂, safe_vector)
    e1 ./= norm(e1)

    e2 = cross(b̂, e1)

    return b̂, e1, e2, B_mag
end

"""
    dipole_field(x, y, z)
    dipole_field(r)
    dipole_field(L, λ)

Compute the dipole magnetic field vector at given position.

Returns the magnetic field components `[Bx, By, Bz]` in tesla, modeled as a magnetic
dipole field. The function can take either three Cartesian coordinates, an array of
Cartesian coordinates, or two components specifying the L-shell and  magnetic latitude of
the position. The field assumes axisymmetry (X-Z plane) and is only valid within 10 Earth
radii from the center of the Earth.

# Arguments

- `x`: Position in x-direction [m].
- `y`: Position in y-direction [m].
- `z`: Position in z-direction [m].
- `r`: Three-element vector `[x, y, z]` of Cartesian coordinates [m].
- `L`: L-shell of the position [RE].
- `λ`: Magnetic latitude of the position [rad].

# Returns

- `SVector{3, Float64}`: Magnetic field components `(Bx, By, Bz)` [T].

# Throws

- `ArgumentError`: If the position is at the origin.
- `ArgumentError`: If the position is inside or on the Earth's surface (`r ≤ RE`).
- `ArgumentError`: If the position is outside the valid range (`r > 10 RE`).
- `ArgumentError`: If the input vector `r` does not have exactly three elements.
"""
function dipole_field(x::Real, y::Real, z::Real)

    r2 = x^2 + y^2 + z^2       #[m]
    r = sqrt(r2)

    iszero(r) && throw(ArgumentError("Dipole field not defined in position origo"))

    r ≤ RE && throw(
        ArgumentError(
            "The position is inside/on the Earth's surface, is this as intended?"
        )
    )

    r > 10 * RE && throw(
        ArgumentError(
            "Position outside of valid range for dipole field approximation, is this as" *
            "intended?"
        )
    )

    r5 = r2^2 * r

    C = - (mu0 / (4π)) * M

    Bx = C * ((3 * x * z) / r5)
    By = C * ((3 * y * z) / r5)
    Bz = C * (3 * (z^2) - r2) / (r5)

    return SVector(Bx, By, Bz)
end

function dipole_field(r::AbstractVector)
    length(r) == 3 || throw(
        ArgumentError("Position vector needs to have three cartesian components, whats up?")
    )

    return dipole_field(r...)
end

function dipole_field(pos::Cartesian)
    return dipole_field(pos.x, pos.y, pos.z)
end

function dipole_field(pos::Dipolar)
    r = Cartesian(pos)
    return dipole_field(r.x, r.y, r.z)
end

"""
    tsyganenko_field(
    x::Real, y::Real, z::Real;
    t="2020-01-01T02:00:00",
    pdyn=2.0,
    dst=-28.0,
    byimf=1.0,
    bzimf=2.0,
    ps=`nothing`,
    component="both",
)
    tsyganenko_field(r::AbstractVector; kwargs...)
    tsyganenko_field(position::Dipolar; kwargs...)
    tsyganenko_field(position::Cartesian; kwargs...)

# Arguments
- `x, y, z`: Cartesian position components [m].
- `r`: Cartesian position vector [m].
- `position`: Position as `Cartesian(x, y, z)` [m] or `Dipolar(L, λ, ϕ)` [RE, rad, rad].

# Keyword Arguments

- `t`: Given as a string in ISO8601 format, default is '"2020-01-01T00:01:40"'.
- `pdyn`: Solar wind dynamic pressure [nPa], default is `2.0`.
- `dst`: Disturbance Storm Time index [nT], default is `-50.0`.
- `byimf`: IMF By component [nT], default is `1.0`.
- `bzimf`: IMF Bz component [nT], default is `2.0`.
- `ps`: Dipole tilt angle, uses `t` to calculate it if the value is `nothing`. The
  default is `nothing`.
- `component`: String to chose which component of the field to return; either external
  `ext`, internal `int` or the total field `both`, default is `"both"`.

# Returns

- A three-element vector `[Bx, By, Bz]` containing magnetic-field components [T].

# Throws

- `ArgumentError`: If invalid argument is given to `component`.
"""
function tsyganenko_field_gsm(
    x::Real,
    y::Real,
    z::Real,
    t::DateTime,
    ps::Real;
    pdyn::Real,
    dst::Real,
    byimf::Real,
    bzimf::Real,
    component::String,
)

    r_RE = [x/RE, y/RE, z/RE]
    param = (; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

    if component == "both"
        B_ext_nT = TS04(param)(r_RE, ps)
        B_int_nT = TsyIGRF()(r_RE, t)
        B_tot_nT = B_ext_nT .+ B_int_nT     # [nT]
        Bx, By, Bz = B_tot_nT .* 1e-9       # [T]
    elseif component == "ext"
        B_ext_nT = TS04(param)(r_RE, ps)  # [nT]
        Bx, By, Bz = B_ext_nT .* 1e-9          # [T]
    elseif component == "int"
        B_int_nT = TsyIGRF()(r_RE, t)       # [nT]
        Bx, By, Bz = B_int_nT .* 1e-9       # [T]
    else
        throw(ArgumentError("Invalid `component` chosen, must be either `ext`, `int` or
        `both`."))
    end

    return SVector(Bx, By, Bz)
end

# Helper function for rotation
function rotate_y(v, α)
    return [cos(α) * v[1] + sin(α) * v[3], v[2], -sin(α) * v[1] + cos(α) * v[3]]
end

function (model::TsyganenkoMagneticField)(r::Cartesian)
    r_mag = sqrt(r.x^2 + r.y^2 + r.z^2)
    r_mag > model.r_limit && error(
        "Tsyganenko evaluated past the limit, r = $(round(r_mag / RE, digits=2)) RE"
    )

    # SM → GSM
    r_gsm = rotate_y(SVector(r.x, r.y, r.z), model.ps)

    # Evaluate field
    B_gsm = tsyganenko_field_gsm(
        r_gsm...,
        model.t,
        model.ps;
        pdyn = model.pdyn,
        dst = model.dst,
        byimf = model.byimf,
        bzimf = model.bzimf,
        component = model.component,
    )

    # GSM → SM
    return rotate_y(B_gsm, -model.ps)
end
