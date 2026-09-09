using AURORA
using LinearAlgebra: norm, cross

"""
    boris_mover_TOF(
    magnetic_field,
    r0,
    v0,
    r_source;
    n_T::Int=100_000,
    resolution::Int=10,
    store_trajectory::Bool=false
)

Performs a numerical calculation on electrons propagating in a magnetic field.

Using a boris-mover, the propagation of an electron under the influence `magnetic field`, a
defined magnetic field function is numerically calculated, from the initial
position `r0` (radial distance from the center of the Earth), with some initial velocity
`v0`, until reaching the desired position `r_source`. The method propagates the particle
away from the Earth to some source, as the available information about the electron is given
at the top of the defined ionosphere.

# Arguments

- `magnetic_field`: Magnetic field function `f(x, y, z)`.
- `r0`: Initial position of the electron [m].
- `v0`: Initial velocity of the electron [m/s].
- `r_source`: The source used for the electron transport, i.e. the stopping condititon for
  the simulation, given as radial distance from the center of the Earth [m].

# Keyword Arguments

- `n_T`: Total number of gyroperiods simulated for. Default is 100_000.
- `resolution`: Number of samples per gyroperiod used, adjusted along with the magnetic
  field strength. Default is 10.
- `store_trajectory`: Option to store the full trajectory of the particle. This requires a
  lot of memory, and is only used on small datasets and for debugging. Default is `false`.

# Returns

- The time of flight, `tof`, for the particle and the footpoint, `(x, y, z)` by default.
  Also returns the full trajectory `r` if `store_trajectory` is `true`.

# Throws

- `ArgumentError`: If `n_T` is zero or negative.
- `ArgumentError`: If `resolution` is zero or negative.

"""
function boris_mover_TOF(
    magnetic_field,
    r0,
    v0,
    r_source;
    n_T::Int=100_000,
    resolution::Int=10,
    store_trajectory::Bool=false
)

    n_T > 0 || throw(ArgumentError("n_T must be positive"))
    resolution > 0 || throw(ArgumentError("Resolution must be positive"))

    # Find total number of steps reqired
    steps = Int(n_T * resolution)

    # Initial position, velocity and TOF
    x, y, z = r0
    vx, vy, vz = v0
    tof = 0.0

    if store_trajectory
        r = zeros(steps+1, 3)
        r[1, :] .= (x, y, z)
    end

    # Update particle
    for i in 2:(steps + 1)
        Bx, By, Bz = magnetic_field(x, y, z)

        # Find desired resolution for current step
        B_mag = sqrt(Bx^2 + By^2 + Bz^2)
        ω_g = gyro_frequency(B_mag, qₑ, mₑ)
        T_g = 2π / ω_g

        # Ensuring no steps are too large
        dt_gyro = abs(T_g / resolution)
        dt_radial = 0.0001 * RE / sqrt(vx^2 + vy^2 + vz^2)
        dt = min(dt_gyro, dt_radial)

        # Half electric acceleration
        q_prime = dt * qₑ / (2mₑ)

        #-----------Core Boris-scheme-----------

        # Magnetic rotation vector
        hx = q_prime * Bx
        hy = q_prime * By
        hz = q_prime * Bz

        h_mag2 = hx^2 + hy^2 + hz^2

        sx = 2hx / (1 + h_mag2)
        sy = 2hy / (1 + h_mag2)
        sz = 2hz / (1 + h_mag2)

        # Boris rotation
        ux_prime = vx + (vy * hz - vz * hy)
        uy_prime = vy + (vz * hx - vx * hz)
        uz_prime = vz + (vx * hy - vy * hx)

        vx += uy_prime * sz - uz_prime * sy
        vy += uz_prime * sx - ux_prime * sz
        vz += ux_prime * sy - uy_prime * sx

        # Update position
        x += vx * dt
        y += vy * dt
        z += vz * dt

        tof += abs(dt)

        if store_trajectory
            r[i, :] .= x, y, z
        end

        # TODO: Add correct check of r_source here
        if z < 1e3
            result = (tof=tof, footpoint = (x, y, z))
            return store_trajectory ? (result..., r=r[1:i,:]) : result
        end
    end

    result = (tof=tof, footpoint=(x, y, z))
    return store_trajectory ? (result..., r=r) : result
end
