using AURORA
using Geodesy
using PythonCall
using Dates

## Coordinate types
struct Cartesian
    x::Float64
    y::Float64
    z::Float64
end

# NOTE: This might be an issue, as it is spherical dipole coordinates, not pure spherical coordinates!
struct Dipolar
    L::Float64
    λ::Float64
    ϕ::Float64
end

struct EastNorthUp
    e::Float64
    n::Float64
    u::Float64
end

struct LTPReference
    lat_deg::Float64
    lon_deg::Float64
    alt_m::Float64
end

struct LatitudeLongitudeAltitude
    lat_deg::Float64
    lon_deg::Float64
    alt_m::Float64
end


# Adding some locations that can be called by name
LTPReference(name::Symbol) = LTPReference(
    Dict(
        :tromso       => (69.66,  18.94, 0.0),
        :longyearbyen => (78.22,  15.63, 0.0),
        :sodankyla    => (67.37,  26.63, 0.0),
    )[name]...
)

# Define conversion between the types
function Dipolar(r::Cartesian)
    r_mag = sqrt(r.x^2 + r.y^2 + r.z^2)
    λ     = asin(r.z / r_mag)
    ϕ     = mod2pi(atan(r.y, r.x))
    L     = r_mag / (RE * cos(λ)^2)
    return Dipolar(L, λ, ϕ)
end

function Cartesian(r::Dipolar)
    r_mag = r.L * RE * cos(r.λ)^2
    x     = r_mag * cos(r.λ) * cos(r.ϕ)
    y     = r_mag * cos(r.λ) * sin(r.ϕ)
    z     = r_mag * sin(r.λ)
    return Cartesian(x, y, z)
end

#===

MAKING COORDINATE TRANSFORMS

===#
# Solar coefficients and subsolar point is found using code from
# Laundal and Richmond (2017), appendix C
function solar_coefficients(year::Int)
    yr = year - 2000
    leapdays = floor(Int, (year - 1601) / 4) - 99
    if year ≤ 1900
        leapdays += 3 - floor(Int, (year - 1601) / 100)
    end

    l0 = -79.549 + (-0.238699  * (yr - 4 * leapdays) + 3.08514e-2 * leapdays)
    g0 = -2.472  + (-0.2558905 * (yr - 4 * leapdays) - 3.79617e-2 * leapdays)

    return (yr=yr, leapdays=leapdays, l0=l0, g0=g0)
end

function subsolar_point(year::Int, doy::Int, ut::Real)
    # Get the solar coefficients
    c = solar_coefficients(year)


    df = (ut/86400 - 1.5) + doy        # Days since the reference epoch
    l_deg = c.l0 + 0.9856474 * df      # mean longitude of the sun
    g_deg = c.g0 + 0.9856003 * df      # mean anomaly of the sun

    # Ecliptic longitude of the sun
    g_rad = deg2rad(g_deg)
    λ_deg = l_deg + 1.915 * sin(g_rad) + 0.020 * sin(2*g_rad)
    sin_λ = sind(λ_deg)

    # Obliquity of the ecliptic
    n = df + 365 * c.yr + c.leapdays
    ε_deg = 23.439 - 4.0e-7 * n

    # Right ascension and declination
    α = rad2deg(atan(cosd(ε_deg) * sin_λ, cosd(λ_deg)))
    δ = rad2deg(asin(sind(ε_deg) * sin_λ))

    # Subsolar colatitude
    subsolar_colatitude = 90 - δ

    # Equation of time wrapped to (-180, 180]
    equation_of_time_deg = l_deg - α
    equation_of_time_deg -= 360 * round(equation_of_time_deg / 360)

    # Sub-solar longitude from the solar time
    subsolar_longitude = 180 - ((ut / 240) + equation_of_time_deg)
    subsolar_longitude -= 360 * round(subsolar_longitude / 360)

    return (colat = subsolar_colatitude, lon = subsolar_longitude)
end

function subsolar_point(t::DateTime)
    ut = hour(t) * 3600 + minute(t) * 60 + second(t)
    return subsolar_point(year(t), dayofyear(t), ut)
end

function sun_unit_vector_GEO(t::DateTime)
    p = subsolar_point(t)
    return [
        sind(p.colat) * cosd(p.lon),
        sind(p.colat) * sind(p.lon),
        cosd(p.colat),
    ]
end


# TODO: Plot for verification
function axis_SM(t::DateTime)
    m̂ = dipole_axis_GEO()               # m̂, unit vector to the north geomagnetic pole
    ŝ = sun_unit_vector_GEO(t)          # ŝ, connects the sun and the Earth center

    ẑ = m̂
    ŷ = normalize(cross(m̂, ŝ))
    x̂ = cross(ŷ, ẑ)
    return x̂, ŷ, ẑ
end


function dipole_tilt(t::DateTime)
    _, _, ẑ = axis_SM(t)
    ŝ = sun_unit_vector_GEO(t)
    ps = asin(dot(ẑ, ŝ))
    return ps
end

function dipole_axis_GEO(;year::Int=2020)
    # Coefficients from IGRF-14
    if year == 2020
        g10, g11, h11 = -29403.41, -1451.37, 4653.35
    elseif year == 2025
        g10, g11, h11 = -29350.0, -1410.3, 4545.5
    else
        throw(ArgumentError("You have to add IGRF-coefficients to use this year..."))
    end

    m = [g11, h11, g10]
    return -m / norm(m)
end

function LLA_to_SM(lat_deg, lon_deg, alt_m, t::DateTime)

    # Take LLA to ECEF using Geodesy
    ref_lla = LLA(lat=lat_deg, lon=lon_deg, alt=alt_m)
    ecef = ECEF(ref_lla, wgs84)

    # Make into Cartesian coordinates
    ref_geo = [ecef.x, ecef.y, ecef.z]

    # Rescale to simulation sphere
    ref_geo_scaled = (RE + z_ionosphere) .* ref_geo ./ norm(ref_geo)
    #ref_geo_scaled = RE .* ref_geo ./ norm(ref_geo)

    # Project onto the SM axes at the correct time
    x̂, ŷ, ẑ = axis_SM(t)
    ref_sm = [dot(ref_geo_scaled, x̂), dot(ref_geo_scaled, ŷ), dot(ref_geo_scaled, ẑ)]

    return ref_sm
end

function EastNorthUp(r::Cartesian, ref::LTPReference)
    # Make ECEF Geodesy-element of position
    ecef = ECEF(r.x, r.y, r.z)

    # Make ECEF Geodesy-element of reference
    ref_lla = LLA(ref.lat_deg, ref.lon_deg, ref.alt)
    ref_ecef = ECEF(ref_lla, wgs84)

    # Convert position to ENU with reference
    enu = ENU(ecef, ref_ecef, wgs84)

    return EastNorthUp(enu.e, enu.n, enu.u)
end

function EastNorthUp(r::Cartesian, ref::Symbol)
    return EastNorthUp(r, LTPReference(ref))
end

function EastNorthUp(r::Dipolar, ref::LTPReference)
    return EastNorthUp(Cartesian(r), ref)
end

function EastNorthUp(r::Dipolar, ref::Symbol)
    return EastNorthUp(Cartesian(r), LTPReference(ref))
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
