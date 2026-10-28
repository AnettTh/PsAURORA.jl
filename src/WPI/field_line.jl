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



"""
    find_R_max(
    magnetic_field,
    position::Cartesian;
    ds=RE*0.01,
    n_steps=10_000,
    store_trace::Bool=false
)
    find_R_max(magnetic_field::AbstractMagneticField, position::Spherical; kwargs...)

Find maximum radial distance from the Earth along a magnetic field line.

The function takes  either an arbitrary magnetic field and a position, or a `ParticleState`,
and traces the field line from some initial position until reaching the ionosphere, in both '
directions. Then, it returns the largest value, as well as the trace if specified.

# Arguments

- `magnetic_field`: Magnetic field function `f(x, y, z)`.
- `position`: Position of particle either `Cartesian(x, y, z)` or `Spherical(L, λ, ϕ)`.

# Keyword Arguments

- `ds`: Size of the step-length along the magnetic field lines, default is `0.01RE` [m].
- `n_steps`: Maximum number of steps taken before ending the tracing, default is `10 000`.
- `store_trace`: Option to store the total trace, default is `false`.

# Throws

- `ArgumentError`: If the magnetic field gives a zero-value, meaning that the initial
  conditions are invalid.

"""
function find_R_max(
    magnetic_field::AbstractMagneticField,
    position::Cartesian;
    ds=RE*0.01,
    n_steps=10_000,
    store_trace::Bool=false
)

    # Initial postition
    x0 = position.x
    y0 = position.y
    z0 = position.z

    # Initial guess
    R_max = norm([x0, y0, z0])

    # If trace is to be returned, split it in forward and backward tracing
    xs_fwd = store_trace ? [x0] : Float64[]
    ys_fwd = store_trace ? [y0] : Float64[]
    zs_fwd = store_trace ? [z0] : Float64[]
    xs_bwd = store_trace ? [x0] : Float64[]
    ys_bwd = store_trace ? [y0] : Float64[]
    zs_bwd = store_trace ? [z0] : Float64[]

    # Trace along fieldlines in both direcitons
    for sign in [1, -1]
        x, y, z = x0, y0, z0
        for step in 1:n_steps
            # Find the current magnetic field direction
            B = magnetic_field(Cartesian(x, y, z))
            B_mag = norm(B)
            iszero(B_mag) && throw(ArgumentError("Why is your magnetic field zero, Miss??"))
            bx, by, bz = B ./ B_mag

            # Move one step in the right direction
            x += sign * ds * bx
            y += sign * ds * by
            z += sign * ds * bz

            # Update the largest sampled radius
            R_current = sqrt(x^2 + y^2 + z^2)
            R_max = max(R_current, R_max)

            # Put in correct container
            if store_trace
                if sign == 1
                    push!(xs_fwd, x); push!(ys_fwd, y); push!(zs_fwd, z)
                else
                    push!(xs_bwd, x); push!(ys_bwd, y); push!(zs_bwd, z)
                end
            end

            # Stop at the ionosphere or above 10 RE
            #R_current ≤ RE + z_ionosphere && break
            if step > 10
                R_current ≤ RE + z_ionosphere && break
            end
            R_current > 10 * RE && break
        end
    end

    if store_trace
        # Assemble and remove duplicates: ionosphere 1 → equator → ionosphere 2
        xs = [reverse(xs_bwd); xs_fwd[2:end]]   # ← remove first element of fwd (duplicate x0)
        ys = [reverse(ys_bwd); ys_fwd[2:end]]
        zs = [reverse(zs_bwd); zs_fwd[2:end]]

        return R_max, (; xs, ys, zs)
    else
        return R_max
    end
end


function find_R_max(magnetic_field::AbstractMagneticField, position::Spherical; kwargs...)
    return find_R_max(magnetic_field, Cartesian(position); kwargs...)
end

# TODO: Make this make sense
function find_L_shell(
    magnetic_field::AbstractMagneticField,
    lat::Float64,
    lon::Float64;
    ds::Float64=RE*0.005,
    n_steps::Int=100_000,
)
    # Radius of the ionsphere with Earth in center
    r_iono = RE + z_ionosphere

    x0 = r_iono * cos(lat) * cos(lon)
    y0 = r_iono * cos(lat) * sin(lon)
    z0 = r_iono * sin(lat)

    function trace_fieldline(x0, y0, z0)
        xs, ys, zs = [x0], [y0], [z0]
        x, y, z = x0, y0, z0

        # Determine outward direction
        B0 = magnetic_field(Cartesian(x0, y0, z0))
        b̂0 = B0 ./ norm(B0)
        r̂0 = [x0, y0, z0] ./ norm([x0, y0, z0])
        sign = dot(b̂0, r̂0) < 0 ? -1 : 1

        for _ in 1:n_steps
            B     = magnetic_field(Cartesian(x, y, z))
            B_mag = norm(B)
            bx, by, bz = B ./ B_mag

            x += sign * ds * bx
            y += sign * ds * by
            z += sign * ds * bz

            push!(xs, x); push!(ys, y); push!(zs, z)

            sqrt(x^2 + y^2 + z^2) > 10 * RE && break
        end
        return xs, ys, zs
    end

    xs, ys, zs = trace_fieldline(x0, y0, z0)

    # Find magnetic field at every point
    Bs = zeros(length(xs))
    for (i, (x, y, z)) in enumerate(zip(xs, ys, zs))
        Bs[i] = norm(magnetic_field(Cartesian(x, y, z)))
    end


    # Find point with lowest fieldstrength
    idx_eq = argmin(Bs)
    r_eq = sqrt(xs[idx_eq]^2 + ys[idx_eq]^2 + zs[idx_eq]^2)
    L = r_eq / RE

    return L
end
