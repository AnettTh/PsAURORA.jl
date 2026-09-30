using AURORA

# NOTE: This isn't even used anywhere?
"""
    dipolar_r_to_λL(r::AbstractVector)
    dipolar_r_to_λL(x::Real, y::Real, z::Real)

Converts a position in a dipolar field from cartesian coordinates to magnetic latitude and
L-shell.

# Arguments

- `r`: Cartesian position (vector format) [m].
- `x, y, z`: Cartesian position (component format) [m].

# Returns

- `λ`: Position magnetic latitude [rad].
- `L`: Position L-shell.

# Throws

- `ArgumentError`: If `r` don't have exactly three components.
- `ArgumentError`: If `r` isn't in the x-z-plane.
- `ArgumentError`: If `r` has zero magnitude.
- `ArgumentError`: If `r` is inside Earth.
"""
function dipolar_r_to_λL(r::AbstractVector)
    length(r) == 3 || throw(ArgumentError("r must have three components (Chartesian)."))
    iszero(r[2]) || throw(ArgumentError("Assumes x-z-plane, your y is invalid. Check it!"))

    r_mag = norm(r)

    iszero(r_mag) && throw(ArgumentError("r must be nonzero"))

    r_mag ≤ RE && throw(ArgumentError("r is inside Earth"))
    λ = asin(r[3] / r_mag)
    L = r_mag / (RE * cos(λ)^2)
    return λ, L
end

function dipolar_r_to_λL(x::Real, y::Real, z::Real)
    return dipolar_r_to_λL([x, y, z])
end
