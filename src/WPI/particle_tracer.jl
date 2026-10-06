using AURORA
using ProgressMeter

"""
    boris_mover(
    particle::ParticleState;
    n_T::Int=100_000,
    resolution::Int=10,
    store_trajectory::Bool=false
)

Simulate the trajectory of a charged particle using the Boris mover.

Uses the information about the particle given by `ParticleState` to numerically trace its
gyromotion in a magnetic field. # TODO: Add more descriptions here as wave and such is added

# Arguments

- `particle`: Struct holding the state of the particle.

# Keyword Arguments

- `n_T`: Number of (initial) gyroperiods to simulate for, default is `100_000`.
- `resolution`: Steps per gyroperiod, default is `10`.
- `store_trajectory`: Option to store the trajectory, default is `false`.

# Returns

- If the trajectroy is stored, a named tuple containing `position`, `velocity`, `time` and
  `footpoint`, if not, only the time-of-flight `time` and the position of impact with the
  ionsophere, `footpoint` is returned.
"""
function boris_mover(
    particle::ParticleState
    # plasma::PlasmaState,
    # wave::AbstractWave
    ;
    n_T::Int=100_000,
    resolution::Int=10,
    store_trajectory::Bool=false
)

    n_T > 0 || throw(ArgumentError("n_T must be positive"))
    resolution > 0 || throw(ArgumentError("resolution must be positive"))

    # Define where the particle has hit the ionosphere
    r_ionosphere = RE + z_ionosphere

    # Find total number of steps reqired
    steps = Int(n_T * resolution)

    # Find initial E-field, B-field and gyroperiods
    E0 = [0.0, 0.0, 0.0]
    #B0 = particle.magnetic_field(particle.r0...)
    #Ω_e0 = Ω_e[0]
    #T0 = 2π / ω_g0

    if store_trajectory
        r_array = zeros(steps + 1, 3)
        v_array = zeros(steps + 1, 3)
        t_array = zeros(steps + 1)
    end

    # Initial values
    x, y, z = particle.r0
    vx, vy, vz = particle.v0
    t_current = 0.0

    if store_trajectory
        r_array[1, :] .= (x, y, z)
        v_array[1, :] .= (vx, vy, vz)
        t_array[1] = t_current
    end

    # Precompute field interpolator # NOTE: This does not work yet, i think???
    B_interpolated = particle.B_interpolated
    #B_interpolated = isnothing(particle.B_interpolated) ?
    #    get_magnetic_field(particle.magnetic_field, particle.L, particle.ϕ_eq) :
    #    particle.B_interpolated
    s = 0.0

    # Update particle
    for i in 2:(steps + 1)
        # TODO: Is it as simple as adding some function for δB and δE here?
        Bx, By, Bz = B_interpolated(s) #.+ wave.δB
        Ex, Ey, Ez = E0 #+ wave.δE

        B_mag = sqrt(Bx^2 + By^2 + Bz^2)
        b̂x, b̂y, b̂z = Bx/B_mag, By/B_mag, Bz/B_mag

        # Define cyclotron frequency and period for timestep
        Ω_e = (B_mag * abs(qₑ)) / mₑ
        T = 2π / Ω_e

        # Define adaptive timestep
        dt = T / resolution

        # Half electric acceleration
        q_prime = dt * qₑ / 2mₑ

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

        # NOTE: not sure if ds_step should be before or after t_current is updated
        ds_step = abs((vx*b̂x + vy*b̂y + vz*b̂z) * dt)
        s += ds_step

        if store_trajectory
            r_array[i, :] .= (x, y, z)
            v_array[i, :] .= (vx, vy, vz)
            t_array[i] = t_current
        end

        if x^2 + y^2 + z^2 ≤ r_ionosphere^2
            if store_trajectory
                return (
                    position=r_array[1:i, :],
                    velocity=v_array[1:i, :],
                    time=t_array[1:i],
                    footpoint=(x, y, z)
                )
            else
                return (time=t_current, footpoint=(x, y, z))
            end
        end
    end

    if store_trajectory
        return (position=r_array, velocity=v_array, time=t_array, footpoint=nothing)
    else
        return (time=t_current, footpoint=nothing)
    end
end
