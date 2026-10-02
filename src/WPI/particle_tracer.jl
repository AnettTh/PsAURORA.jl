using AURORA

"""

Simulate the trajectory of a charged particle using the Boris mover.

# Arguments

- `magnetic_field`: function that takes position.
- `electric_field`: function that takes position.
- `n_T::Integer`: number of (initial) gyroperiods to simulate for.
- `r0`: initial position.
- `v0`: initial velocity.
- `q`: particle charge.
- `m`: particle mass.

# Keyword Arguments

- `resolution::Integer=10`: steps per gyroperiod.
- `stop_at_ionosphere::Bool=true`: option to stop the simulation at the ionosphere.

# Returns

A named tuple containing `position`, `velocity`, `time` and `footpoint`
"""
function boris_mover(
    magnetic_field::AbstractMagneticField,
    wave::AbstractWave,
    particle::AbstractParticle,
    n_T::Int,
    r0,
    v0,
    q,
    m;
    resolution::Int=10
)

    n_T > 0 || throw(ArgumentError("n_T must be positive"))
    resolution > 0 || throw(ArgumentError("resolution must be positive"))
    m != 0 || throw(ArgumentError("particle mass must be nonzero"))

    # Define where the particle has hit the ionosphere
    r_ionosphere = RE + z_ionosphere
    #footpoint = nothing

    # Find total number of steps reqired
    steps = Int(n_T * resolution)

    # Find initial E-field, B-field and gyroperiod
    #E0 = electric_field(r0...)
    #B0 = magnetic_field(r0...)
    #ω_g0 = gyro_frequency(B0, q, m)
    #T_g0 = 2π / ω_g0

    # Make time-range based on resolution (samples per gyroperiod)
    #dt = T_g0 / resolution

    r_array = zeros(steps + 1, 3)
    v_array = zeros(steps + 1, 3)
    t_array = zeros(steps + 1)

    # Initial values
    x, y, z = r0
    vx, vy, vz = v0
    t_current = 0.0

    r_array[1, :] .= (x, y, z)
    v_array[1, :] .= (vx, vy, vz)
    t_array[1] = t_current

    # Update particle
    for i in 2:(steps + 1)
        Bx, By, Bz = magnetic_field(x, y, z)
        Ex, Ey, Ez = electric_field(x, y, z)

        ω_g = gyro_frequency([Bx, By, Bz], q, m)
        T_g = 2π / ω_g

        # Define adaptive timestep
        dt = T_g / resolution

        # Half electric acceleration
        q_prime = dt * q / 2m

        #-----------Core Boris-scheme-----------

        # Half electric acceleration
        vx_minus = vx + q_prime * Ex
        vy_minus = vy + q_prime * Ey
        vz_minus = vz + q_prime * Ez

        # Magnetic rotation vector
        hx = q_prime * Bx
        hy = q_prime * By
        hz = q_prime * Bz
        h_mag2 = hx^2 + hy^2 + hz^2

        sx = 2hx / (1 + h_mag2)
        sy = 2hy / (1 + h_mag2)
        sz = 2hz / (1 + h_mag2)

        ux_prime = vx_minus + (vy_minus * hz - vz_minus * hy)
        uy_prime = vy_minus + (vz_minus * hx - vx_minus * hz)
        uz_prime = vz_minus + (vx_minus * hy - vy_minus * hx)

        vx_plus = vx_minus + (uy_prime * sz - uz_prime * sy)
        vy_plus = vy_minus + (uz_prime * sx - ux_prime * sz)
        vz_plus = vz_minus + (ux_prime * sy - uy_prime * sx)

        vx = vx_plus + q_prime * Ex
        vy = vy_plus + q_prime * Ey
        vz = vz_plus + q_prime * Ez

        x += vx * dt
        y += vy * dt
        z += vz * dt

        # Update stored values with current values
        t_current += dt

        r_array[i, :] .= (x, y, z)
        v_array[i, :] .= (vx, vy, vz)
        t_array[i] = t_current

        if stop_at_ionosphere && x^2 + y^2 + z^2 ≤ r_ionosphere^2

            return (
                position=r_array[1:i, :],
                velocity=v_array[1:i, :],
                time=t_array[1:i],
                footpoint=(x, y, z)
            )
        end
    end

    return (position = r_array, velocity = v_array, time = t_array, footpoint = nothing)
end
