using AURORA
using LinearAlgebra
using Dates

abstract type AbstractFieldLine end

struct StationFieldLine{F<:AbstractMagneticField} <: AbstractFieldLine
    magnetic_field::F
    # Geometry, all in SM coordinates and [m] or [rad], from south FP → north FP
    xs::Vector{Float64}     # x-positions along the line
    ys::Vector{Float64}     # y-positions along the line
    zs::Vector{Float64}     # z-positions along the line
    rs::Vector{Float64}     # Radius at every point
    λs::Vector{Float64}     # Latitude at every point
    ϕs::Vector{Float64}     # Longitude at every point
    s::Vector{Float64}      # Arc length from equator, negative is southern hemisphere

    Bs::Vector{Float64}     # Fieldstrength |B| at each point along the line [T]
    # Magnetic equator ⟹ point of minimum |B|
    i_eq::Int               # Index of equator
    position_eq::Cartesian  # Cartesian position at equator
    r_eq::Float64           # Radius at equator [m]
    λ_eq::Float64           # Latitude at equator [rad]
    ϕ_eq::Float64           # Longitude at equator [rad]
    B_eq::Float64           # Fieldstrength at equator [T]

    r_max::Float64      # Maximal radius, not necessary the same as equator [m].

    position_fp_south::Cartesian    # Cartesian southern footpoint [m]
    position_fp_north::Cartesian    # Cartesian northern footpoint [m]
    B_along_index::Bool             # Motion of particle relative to indexation.
end


function trace_fieldline(
    magnetic_field::AbstractMagneticField,
    position::Cartesian;
    ds::Float64=0.0005RE,
    n_steps::Int=100_000,
    r_limit::Float64=15RE
)
    # Radius of the ionosphere and initial position
    r_iono = RE + z_ionosphere
    x, y, z = position.x, position.y, position.z

    # Determine outward direction
    B0  = magnetic_field(position)
    dir = (B0[1]*x + B0[2]*y + B0[3]*z) < 0 ? -1 : 1

    # Store values to return
    xs, ys, zs, Bs = [x], [y], [z], Float64[]
    closed_fieldline = false

    # Trace the fieldline
    for _ in 1:n_steps

        # Check the magnetic field
        B  = magnetic_field(Cartesian(x, y, z))
        Bm = norm(B)
        (iszero(Bm) || isnan(Bm)) && throw(
            ArgumentError("Field is zero or NaN at ($x, $y, $z)")
        )

        # Store the current |B| before leaving it
        push!(Bs, Bm)

        # Update position
        x += dir * ds * B[1] / Bm
        y += dir * ds * B[2] / Bm
        z += dir * ds * B[3] / Bm
        push!(xs, x); push!(ys, y); push!(zs, z)

        # Check if the current position is past the limit
        r = sqrt(x^2 + y^2 + z^2)
        r > r_limit && error("Field line did not close within $(r_limit / RE) RE")

        # Make a note if the fieldline is closed and end the trace
        if r ≤ r_iono
            closed_fieldline = true
            break
        end
    end

    # Warn if too few steps are used before reaching the limit or the ionosphere
    closed_fieldline || error(
        "Field line did not return to the ionosphere in $n_steps steps"
    )

    # The last point is just below r_iono, move it onto the sphere and get its |B|
    r_in, r_out = sqrt(xs[end-1]^2 + ys[end-1]^2 + zs[end-1]^2), sqrt(x^2 + y^2 + z^2)
    f = (r_in - r_iono) / (r_in - r_out)
    xs[end] = xs[end-1] + f * (xs[end] - xs[end-1])
    ys[end] = ys[end-1] + f * (ys[end] - ys[end-1])
    zs[end] = zs[end-1] + f * (zs[end] - zs[end-1])
    push!(Bs, norm(magnetic_field(Cartesian(xs[end], ys[end], zs[end]))))

    # Change order such that it spans from south to north
    if zs[end] < zs[1]
        reverse!(xs); reverse!(ys); reverse!(zs); reverse!(Bs)
    end
    return xs, ys, zs, Bs
end


function StationFieldLine(
    magnetic_field,
    ref::LTPReference,
    t::DateTime;
    ds=RE*0.005,
    n_steps=10_000
)
    # Get the initial condition based on reference position and time
    position_sm = LLA_to_SM(ref.lat_deg, ref.lon_deg, 0.0, t)
    position = Cartesian(position_sm...)

    # Trace the fieldline
    xs, ys, zs, Bs = trace_fieldline(magnetic_field, position; ds=ds, n_steps=n_steps)
    rs = @. sqrt(xs^2 + ys^2 + zs^2)
    λs = @. asin(zs ./ rs)
    ϕs = @. mod2pi.(atan.(ys, xs))

    # Get index of equator
    i_eq = argmin(Bs)

    # Direction of particles movement compared to the indexation of the fieldline
    # NOTE: This might not be necessary, check it!
    j = clamp(i_eq, 2, length(xs) - 1)
    tangent = SVector(xs[j+1] - xs[j-1], ys[j+1] - ys[j-1], zs[j+1] - zs[j-1])
    B_along_index = dot(magnetic_field(Cartesian(xs[j], ys[j], zs[j])), tangent) > 0

    # Get arc-length from magnetic equator, negative is towards the south
    total_s = [0.0; cumsum(sqrt.(diff(xs).^2 .+ diff(ys).^2 .+ diff(zs).^2))]
    s = total_s .- total_s[i_eq]

    # Get equatorial values also, using i_eq
    position_eq = Cartesian(xs[i_eq], ys[i_eq], zs[i_eq])
    r_eq = rs[i_eq]
    λ_eq = λs[i_eq]
    ϕ_eq = ϕs[i_eq]
    B_eq = Bs[i_eq]

    # Get r_max
    i_r_max = argmax(rs)
    r_max = rs[i_r_max]

    # Get footpoint positions
    position_fp_south = Cartesian(xs[1], ys[1], zs[1])
    position_fp_north = Cartesian(xs[end], ys[end], zs[end])

    return StationFieldLine(
        magnetic_field,

        xs,
        ys,
        zs,
        rs,
        λs,
        ϕs,
        s,

        Bs,

        i_eq,
        position_eq,
        r_eq,
        λ_eq,
        ϕ_eq,
        B_eq,

        r_max,

        position_fp_south,
        position_fp_north,
        B_along_index
        )
end

function Base.show(io::IO, line::StationFieldLine)
    # Helpers that work on the Cartesian footpoints, so they don't depend on the array fields
    radius(p)   = sqrt(p.x^2 + p.y^2 + p.z^2) / RE
    latitude(p) = rad2deg(asin(p.z / sqrt(p.x^2 + p.y^2 + p.z^2)))
    mlt(p)      = longitude_to_MLT(atan(p.y, p.x); degrees=false)
    fmt(x, d)   = round(x, digits=d)

    fp_s, fp_n = line.position_fp_south, line.position_fp_north
    length_RE  = (line.s[end] - line.s[1]) / RE

    lines = [
        "StationFieldLine ($(nameof(typeof(line.magnetic_field))))",
        "    Trace length: $(fmt(length_RE, 4)) RE",
        "    r_eq:         $(fmt(line.r_eq / RE, 4)) RE",
        "    λ_eq:         $(fmt(rad2deg(line.λ_eq), 4))°",
        "    ϕ_eq:         $(fmt(rad2deg(line.ϕ_eq), 4))° ($(fmt(longitude_to_MLT(line.ϕ_eq; degrees=false), 2)) MLT)",
        "    B_eq:         $(fmt(line.B_eq * 1e9, 4)) nT",
        "    r_max:        $(fmt(line.r_max / RE, 4)) RE",
        "    FP south:     λ = $(fmt(latitude(fp_s), 4))°, $(fmt(mlt(fp_s), 4)) MLT, r = $(fmt(radius(fp_s), 3)) RE",
        "    FP north:     λ = $(fmt(latitude(fp_n), 4))°, $(fmt(mlt(fp_n), 4)) MLT, r = $(fmt(radius(fp_n), 3)) RE",
    ]
    print(io, join(lines, "\n"))
end
