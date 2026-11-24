using AURORA
using Geodesy
using PythonCall
using Dates

#-------------------------------------Coordinate types--------------------------------------
"""
    Cartesian

Cartesian coordinates, given in metres. In the code, it is used in SM coordinates, unless
stated otherwise.

# Fields
# TODO: Check that this is correct!
- `x`: Positive towards the sun [m].
- `y`: Positive towards dusk [m].
- `z`: Positive north along the dipole axis [m].

# Constructors

    Cartesian(x::Float64, y::Float64, z::Float64)
"""
struct Cartesian
    x::Float64
    y::Float64
    z::Float64
end


"""
    EastNorthUp

Cartesian coordinates centered at some reference point with a local tangent plane, pointing
east, north and up, given in metres.

# Fields

- `e`: Component pointing eastward from the reference [m].
- `n`: Component pointing northward from the reference [m].
- `u`: Component pointing upward from the reference [m].
"""
struct EastNorthUp
    e::Float64
    n::Float64
    u::Float64
end

"""
    LTPReference

Stores latitude, longitude and altitude of a reference point, used to construct a local
tangent plane, used as the reference for the east-north-up coordinate system. Here, the
coordinates are given in degrees rather than radians, as it is more natural to link them to
actual geographic locations.

# Fields

- `lat_deg`: latitude of the reference-position [°].
- `lon_deg`: longitude of the reference-position [°].
- `alt_m`: Altitude of the reference position [m].

# Constructors

    LTPReference(lat_deg, lon_deg, alt_m)
    LTPReference(name::Symbol)

The form using symbols looks up a station existing in the code, all at altitude 0:
`:longyearbyen`, `:tromso` or `:sodankyla`, throwing a `KeyError` for any name not yet added
to the list.
"""
struct LTPReference
    lat_deg::Float64
    lon_deg::Float64
    alt_m::Float64
end

LTPReference(name::Symbol) = LTPReference(
    Dict(
        :longyearbyen => (78.22, 15.63, 0.0),
        :tromso       => (69.66, 18.94, 0.0),
        :sodankyla    => (67.37, 26.63, 0.0),
    )[name]...
)


#-------------------------------MAKING COORDINATE TRANSFORMS--------------------------------

# Solar coefficients and subsolar point is found using code from Laundal and Richmond (2017)
# , appendix C
"""
    solar_coefficients(year::Int)

Calculates the solar coefficients used further to make the conversions from
latitude-lingitude-altitude coordinates to solar-magnetic coordinates.

# Fields

- `year`: Integer representation of the year to use for the calculation.

# Returns

A named tuple `(yr, leapdays, l0, g0)`:

- `yr`: Years since 2000.
- `leapdays`: Number of leap days counted from 1601, offset so that it is 5 for 2020.
- `l0`: Mean longitude of the Sun at the reference epoch of that year [°].
- `g0`: Mean anomaly of the Sun at the reference epoch of that year [°].

# References

- Laundal and Richmond (2017), appendix C.
"""
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

"""
    subsolar_point(year::Int, doy::Int, ut::Real)

Calculates the subsolar colatitude and subsolar longitude, based on the given year, DOY and
time.

# Arguments

- `year`: Current year as integer.
- `doy`: Current day-of-the-year as an integer.
- `ut`:
"""
function subsolar_point(year::Int, doy::Int, ut::Real)
    # Get the solar coefficients
    c = solar_coefficients(year)

    df = (ut/86400 - 1.5) + doy        # Days since the reference epoch
    l_deg = c.l0 + 0.9856474 * df      # Mean longitude of the sun
    g_deg = c.g0 + 0.9856003 * df      # Mean anomaly of the sun

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
    scl = sind(p.colat)
    return [
        scl * cosd(p.lon),
        scl * sind(p.lon),
        cosd(p.colat),
    ]
end

function axis_SM(t::DateTime)
    m̂ = dipole_axis_GEO()               # m̂, unit vector to the north geomagnetic pole
    ŝ = sun_unit_vector_GEO(t)          # ŝ, connects the sun and the Earth center

    ŷ = normalize(cross(m̂, ŝ))
    x̂ = cross(ŷ, m̂)
    return x̂, ŷ, m̂      # ẑ = m̂
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
        throw(ArgumentError("You have to look up more IGRF-coefficients to use this year"))
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

function SM_to_GEO(r::Cartesian, t::DateTime)
    x̂, ŷ, ẑ = axis_SM(t)
    geo = r.x .* x̂ .+ r.y .* ŷ .+ r.z .* ẑ
    return Cartesian(geo...)
end

function SM_to_ENU(r_sm::Cartesian, ref::LTPReference, t::DateTime)
    # SM → GEO
    r_geo = SM_to_GEO(r_sm, t)

    # Make ECEF Geodesy-element of position
    ecef = ECEF(r_geo.x, r_geo.y, r_geo.z)

    # Make ECEF Geodesy-element of reference
    ref_lla = LLA(ref.lat_deg, ref.lon_deg, ref.alt_m)
    ref_ecef = ECEF(ref_lla, wgs84)

    # Convert position to ENU with reference
    enu = ENU(ecef, ref_ecef, wgs84)

    return EastNorthUp(enu.e, enu.n, enu.u)
end

function SM_to_ENU(r_sm::Cartesian, ref::Symbol, t::DateTime)
    return SM_to_ENU(r_sm, LTPReference(ref), t)
end


function FP_to_ENU(fp, ref::LTPReference, t::DateTime)
    any(isnan, fp) && return [NaN, NaN, NaN]
    enu = SM_to_ENU(Cartesian(fp...), ref, t)
    return [enu.e, enu.n, enu.u]
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
function longitude_to_MLT(ϕ; degrees::Bool=false)
    if !degrees
        return  mod(12.0 + (ϕ / 15.0), 24.0)
    else
        return mod(12.0 + (rad2deg(ϕ) / 15.0), 24.0)
    end
end
