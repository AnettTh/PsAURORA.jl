using AURORA

## Coordinate types
struct Cartesian
    x :: Float64
    y :: Float64
    z :: Float64
end

struct Spherical
    L :: Float64
    λ :: Float64
    ϕ :: Float64
end


# Define conversion between the types
function Spherical(r::Cartesian)
    r_mag = sqrt(r.x^2 + r.y^2 + r.z^2)
    λ     = asin(r.z / r_mag)
    ϕ     = atan(r.y, r.x)
    L     = r_mag / (RE * cos(λ)^2)
    return Spherical(L, λ, ϕ)
end

function Cartesian(r::Spherical)
    r_mag = r.L * RE * cos(r.λ)^2
    x     = r_mag * cos(r.λ) * cos(r.ϕ)
    y     = r_mag * cos(r.λ) * sin(r.ϕ)
    z     = r_mag * sin(r.λ)
    return Cartesian(x, y, z)
end


"""
    longitude_to_MLT(ϕ; degrees::Bool=true)

Convert from longitude to Magnetic Local Time (MLT).

Takes the longitude in either degrees (default) or radians, and converts the position to
magnetic local time.

# Arguments

- `ϕ`: Longitude [° or rad].

# Keyword Arguments

- `degrees`: Indicates the units of the input longitude, default is degrees.
"""
function longitude_to_MLT(ϕ; degrees::Bool=true)
    if degrees
        return  mod(12.0 + (ϕ / 15.0), 24.0)
    else
        return mod(12.0 + (rad2deg(ϕ) / 15.0), 24.0)
    end
end
