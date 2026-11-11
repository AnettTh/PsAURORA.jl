using AURORA
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
