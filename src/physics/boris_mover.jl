using AURORA
using LinearAlgebra: norm, cross

# TODO: Look into making this run faster
function boris_mover_TOF(
    magnetic_field,
    r0,
    v0,
    z_end;
    n_T::Int=100_000,
    resolution::Int=10
)

    n_T > 0 || throw(ArgumentError("n_T must be positive"))
    resolution > 0 || throw(ArgumentError("resolution must be positive"))
    mₑ != 0 || throw(ArgumentError("particle mass must be nonzero"))

    # Find total number of steps reqired
    steps = Int(n_T * resolution)

    # Initial position and velocity
    x, y, z = r0
    vx, vy, vz = v0

    # Find initial B-field and gyroperiod
    # TODO: This might need to be altered when using the specific position
    B0 = magnetic_field(0.0, 0.0, z_end)
    ω_g0 = gyro_frequency(B0, qₑ, mₑ)
    T_g0 = 2π / ω_g0

    # Make time-range based on resolution (samples per gyroperiod)
    dt = abs(T_g0 / resolution)

    # Half electric acceleration
    q_prime = dt * qₑ / (2mₑ)

    # Accumulating time-of-flight
    tof = 0.0
    r = zeros(steps+1, 3)
    r[1, :] .= (x, y, z)

    # Update particle
    for i in 2:(steps + 1)
        Bx, By, Bz = magnetic_field(x, y, z)

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
        r[i, :] .= x, y, z

        if z < 1e3
            return (
                tof,
                footpoint = (x, y, z),
                r,
                i)
        end

    end

    return (
        tof,
        footpoint = (x, y, z),
        r
    )
end
