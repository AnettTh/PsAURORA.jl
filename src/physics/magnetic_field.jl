#using Dates
#using TsyganenkoModels

# TODO: Add tsyganenko-option

"""
    dipole_field(x, y, z)

Make vector containing dipole magnetic field values at given position.

The function takes in Cartesian coordinates (either scalar coordinate or vector) and returns
a vector giving the magnetic field strength corresponding to a modeled dipole field. If the
position is outside of the valid range for the approximation, returns zero-vector for
consistency, with a warning.

# Arguments

- `x`: Position in x-direction, given in meters.
- `y`: Position in y-direction, given in meters.
- `z`: Position in z-direction, given in meters.

# Returns

- A three-element vector `[Bx, By, Bz]` containing magnetic-field components in tesla.

# Throws

- `ArgumentError`: If the position is at the origin or if the vector has the wrong length.
"""
function dipole_field(x, y, z)

    r = sqrt(x^2 + y^2 + z^2)       #[m]

    iszero(r) && throw(ArgumentError("Dipole field not defined in position origo"))

    if r ≤ RE
        @warn "The position is inside/on the Earth's surface, is this as intended?"
    end

    C = - (μ₀ / (4π)) * M

    if r ≤ 7 * RE
        Bx = C * ((3 * x * z) / r^5)
        By = C * ((3 * y * z) / r^5)
        Bz = C * (3 * (z^2) - r^2) / (r^5)
        return [Bx, By, Bz]
    else
        @warn "Position outside of valid range for dipole field approximation, is this as
        intended? Returning zero."

        return [0.0, 0.0, 0.0]
    end
end


function dipole_field(r)
    length(r) == 3 || throw(
        ArgumentError("Position vector needs to have three cartesian components, whats up?")
    )

    return dipole_field(r...)
end


# NOTE: Move elsewhere?
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
