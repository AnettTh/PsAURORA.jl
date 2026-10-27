using AURORA
using StaticArrays

"""
    FieldLine

Precomputed geometry of a magnetic field line, including arc length grid,
Cartesian coordinates, interpolated field, and equatorial reference point.

# Fields
- `B_interpolated`: Interpolated magnetic field as a function of arc length s [T]
- `s_grid`: Arc length grid along the field line [m]
- `s0`: Arc length at the magnetic equator (minimum |B|) [m]
- `r_eq`: Cartesian position at the magnetic equator [m]
- `xs_grid`, `ys_grid`, `zs_grid`: Cartesian coordinates along the field line [m]
"""
struct FieldLine{I}
    B_interpolated :: I
    s_grid         :: Vector{Float64}
    s0             :: Float64
    r_eq           :: SVector{3, Float64}
    B_eq           :: Float64
    xs_grid        :: Vector{Float64}
    ys_grid        :: Vector{Float64}
    zs_grid        :: Vector{Float64}
end

"""
    FieldLine(magnetic_field, L, ϕ_eq)

Calculate trace and do precomputation of magnetic field line at given L-shell and longitude.

Computes the arc length grid, Cartesian coordinates, interpolated field strength and
equatorial reference point (minimum |B| along the field line).

# Arguments
- `magnetic_field`: Magnetic field model.
- `L`: L-shell.
- `ϕ_eq`: Longitude of the field line [rad].

# Returns
"""
function FieldLine(magnetic_field, L, ϕ_eq)

    B_interpolated, s_grid, xs_grid, ys_grid, zs_grid = get_magnetic_field(
        magnetic_field, L, ϕ_eq
    )

    # Find position of lowest magnetic field strength
    i_eq = argmin([norm(B_interpolated(s)) for s in s_grid])
    s0   = s_grid[i_eq]
    r_eq = SVector(xs_grid[i_eq], ys_grid[i_eq], zs_grid[i_eq])
    B_eq = norm(magnetic_field(Cartesian(r_eq...)))

    return FieldLine(B_interpolated, s_grid, s0, r_eq, B_eq, xs_grid, ys_grid, zs_grid)
end
