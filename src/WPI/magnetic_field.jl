#using Dates
#using TsyganenkoModels
using StaticArrays

# TODO: Add tsyganenko-option

"""
    dipole_field(x, y, z)
    dipole_field(r)
    dipole_field(L, λ)

Compute the dipole magnetic field vector at given position.

Returns the magnetic field components `[Bx, By, Bz]` in tesla, modeled as a magnetic
dipole field. The function can take either three Cartesian coordinates, an array of
Cartesian coordinates, or two components specifying the L-shell and  magnetic latitude of
the position. The field is only valid within 10 Earth radii from the center of the Earth.

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
function dipole_field(x, y, z)

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
            "Position outside of valid range for dipole field approximation, is this as
            intended?"
        )
    )

    r5 = r2^2 * r

    C = - (μ₀ / (4π)) * M

    Bx = C * ((3 * x * z) / r5)
    By = C * ((3 * y * z) / r5)
    Bz = C * (3 * (z^2) - r2) / (r5)

    return SVector(Bx, By, Bz)
end


function dipole_field(r)
    length(r) == 3 || throw(
        ArgumentError("Position vector needs to have three cartesian components, whats up?")
    )

    return dipole_field(r...)
end

# TODO: Some throws/warnings?
function dipole_field(L::Real, λ::Real)
    r = L * RE * cos(λ)^2
    x = r * cos(λ)
    z = r * sin(λ)

    return dipole_field(x, 0.0, z)
end

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


# TODO: Doc-string, throws, multiple dispatch
"""
    r_to_λL(r)

Converts a position in a dipolar field from cartesian coordinates to magnetic latitude and
L-shell.

# Arguments

- `r`: Cartesian position [m].

# Returns

- `λ`: Position magnetic latitude [rad].
- `L`: Position L-shell.

# Throws

- `ArgumentError`: If `r` don't have exactly three components.
- `ArgumentError`: If `r` isn't in the x-z-plane. # TODO: Fix this??
- `ArgumentError`: If `r` has zero magnitude.
- `ArgumentError`: If `r` is inside Earth.
"""
function r_to_λL(r)
    length(r) == 3 || throw(ArgumentError("r must have three components (Chartesian)."))
    iszero(r[2]) || throw(ArgumentError("Assumes x-z-plane, your y is invalid. Check it!"))

    r_mag = norm(r)

    iszero(r_mag) && throw(ArgumentError("r must be nonzero"))

    r_mag ≤ RE && throw(ArgumentError("r is inside Earth"))
    λ = asin(r[3] / r_mag)
    L = r_mag / (RE * cos(λ)^2)
    return λ, L
end


#"""
#    tsyganenko_field(
#    x,
#    y,
#    z;
#    time="2020-01-01T00:01:40",
#    pdyn=2.0,
#    dst=-87.0,
#    byimf=2.0,
#    bzimf=-5.0,
#)
#
## Arguments
#- `x`: Position in x-direction, given in meters.
#- `y`: Position in y-direction, given in meters.
#- `z`: Position in z-direction, given in meters.
#
## Keyword Arguments
#
#- `time`: (default "2020-01-01T00:01:40").
#- `pdyn`: Solar wind dynamic pressure [nPa] (default 2.0).
#- `dst`: Disturbance Storm Time index (default -87.0).
#- `byimf`: IMF By component (default 2.0).
#- `bzimf`: IMF Bz component (default -5.0).
#- `component`: Chose which component of the field to return; either external `ext`, internal
#               `int` or the total field `both`.
#
## Returns
#
#- A three-element vector `[Bx, By, Bz]` containing magnetic-field components in tesla.
#
## Throws
#
#- `ArgumentError`: If invalid argument is given to `component`.
#"""
#function tsyganenko_field(
#    x,
#    y,
#    z;
#    time="2020-01-01T00:01:40",
#    pdyn=2.0,
#    dst=-87.0,
#    byimf=2.0,
#    bzimf=-5.0,
#    component="both",
#)
#    t = DateTime(time)
#    r_RE = [x/RE, y/RE, z/RE]
#
#    param = (; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)
#
#    # Construct the external- and internal magnetic field, which then is combined
#    B_ext_nT = TS04(param)(r_RE, ps)
#    B_int_nT = TsyIGRF()(r_RE, t)
#
#
#    if component == "both"
#        B_tot_nT = B_ext_nT .+ B_int_nT
#        Bx, By, Bz = B_tot_nT .* 1e-9
#    elseif component == "ext"
#        Bx, By, Bz = B_ext_nT .* 1e-9
#    elseif component == "int"
#        Bx, By, Bz = B_int_nT .* 1e-9
#    else
#        throw(ArgumentError("Invalid `component` chosen, must be either `ext`, `int` or
#        `both`."))
#    end
#
#
#    return [Bx, By, Bz]
#end
#
