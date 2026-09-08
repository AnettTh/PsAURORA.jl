using AURORA
using LinearAlgebra: norm, cross

# TODO: Look into making this run faster
function boris_mover_TOF(
    magnetic_field,
    r0,
    v0,
    z_end;
    n_T::Int=100_000,
    resolution::Int=10,
    store_trajectory::Bool=false
)

    n_T > 0 || throw(ArgumentError("n_T must be positive"))
    resolution > 0 || throw(ArgumentError("resolution must be positive"))
    mₑ != 0 || throw(ArgumentError("particle mass must be nonzero"))

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
        # TODO: Add multiple-dispatch to use B_mag here
        B_mag = Bx^2 + By^2 + Bz^2
        ω_g = gyro_frequency(B_mag, qₑ, mₑ)
        T_g = 2π / ω_g

        # Ensuring no steps are too large
        dt_gyro = abs(T_g / resolution)
        dt_radial = 0.01 * RE / sqrt(vx^2 + vy^2 + vz^2)
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

        # TODO: Add correct check of z_end here
        if z < 1e3
            result = (tof=tof, footpoint = (x, y, z), steps_taken=i)
            return store_trajectory ? (result..., r=r[1:i,:]) : result
        end
    end

    result = (tof=tof, footpoint=(x, y, z))
    return store_trajectory ? (result..., r=r) : result
end
