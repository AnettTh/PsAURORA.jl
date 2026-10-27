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
gyromotion in a magnetic field.

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
    particle::ParticleState;
    wave::Union{WaveState, Nothing}=nothing,
    λ_start=nothing,
    φ0=0.0,
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
    B0 = particle.B0
    #Ω_e0 = Ω_e[1]
    #T0 = 2π / ω_g0

    if store_trajectory
        r_array = zeros(steps + 1, 3)
        v_array = zeros(steps + 1, 3)
        t_array = zeros(steps + 1)
        α_array = zeros(steps + 1)
    end

    # Figure out where to start the particle from
    if isnothing(λ_start)
        # Use default positions, usually in the equatorial plane
        x, y, z = particle.r0
        vx, vy, vz = particle.v0
        s = particle.field_line.s0
    else
        # Get position of resonance
        r_start = Cartesian(Spherical(particle.L, λ_start, particle.ϕ_eq))
        x, y, z = r_start.x, r_start.y, r_start.z

        # Find v0 at resonance latitude: Evaluate the exact magnetic field once
        B_start = particle.magnetic_field(Spherical(particle.L, λ_start, particle.ϕ_eq))
        b̂, e1, e2, B_start_mag = magnetic_basis(B_start)

        α_start = pitch_angle_at_λ(particle.α_eq, particle.B_eq, B_start_mag)
        v_par   = - particle.v * cos(α_start) * sign(particle.μ)
        v_perp  = particle.v * sin(α_start)

        vx, vy, vz = v_par .* b̂ .+ v_perp * cos(φ0) .* e1 .+ v_perp * sin(φ0) .* e2

        # Get first value of s to represent the start latitude
        # NOTE: If-else-logic can be removed once this tests successful
        if !isnothing(particle.field_line)
            # Find the point in s closest to the initial position
            dists2 = (
                (particle.field_line.xs_grid .- x).^2 .+
                (particle.field_line.ys_grid .- y).^2 .+
                (particle.field_line.zs_grid .- z).^2
            )
            s = particle.field_line.s_grid[argmin(dists2)]
        end
    end

    t_current = 0.0

    if store_trajectory
        r_array[1, :] .= (x, y, z)
        v_array[1, :] .= (vx, vy, vz)
        t_array[1] = t_current
        α_array[1] = acos(abs(particle.μ))
    end

    # Initialize phase if wave is to be considered
    Ψ = isnothing(wave) ? 0.0 : wave.Ψ0

    # Update particle
    @showprogress for i in 2:(steps + 1)

        # Background fields
        if isnothing(particle.field_line)
            B = particle.magnetic_field(Cartesian(x, y, z))
            Bx, By, Bz = B[1], B[2], B[3]
        else
            Bx, By, Bz = particle.field_line.B_interpolated(s)
        end

        Ex, Ey, Ez = E0

        # Add wave perturbation at current time-step
        if !isnothing(wave)
            r_mag = sqrt(x^2 + y^2 + z^2)
            λ_current = asin(z / r_mag)

            # Compute current s
            if !isnothing(particle.field_line)
                dists2 = (particle.field_line.xs_grid .- x).^2 .+
                         (particle.field_line.ys_grid .- y).^2 .+
                         (particle.field_line.zs_grid .- z).^2
                s_current = particle.field_line.s_grid[argmin(dists2)]
            else
                s_current = 0.0   # placeholder — only matters for amplitude models that depend on s
            end

            δB, δE = wave_fields(wave, t_current, s_current, Ψ, SVector(Bx, By, Bz), λ_current)

            Bx += δB[1]
            By += δB[2]
            Bz += δB[3]

            Ex += δE[1]
            Ey += δE[2]
            Ez += δE[3]
        end

        B_mag = sqrt(Bx^2 + By^2 + Bz^2)
        b̂x, b̂y, b̂z = Bx/B_mag, By/B_mag, Bz/B_mag

        # Define cyclotron frequency and period for timestep
        Ω_e = (B_mag * abs(qₑ)) / mₑ
        T = 2π / Ω_e

        # Define adaptive timestep
        dt = T / resolution

        # Half electric acceleration
        q_prime = dt * (-qₑ) / 2mₑ      # NOTE: qₑ here is signed, should this be done in the constant-file instead? Check consistency!

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

        # Update velocity
        vx = vx_plus + q_prime * Ex
        vy = vy_plus + q_prime * Ey
        vz = vz_plus + q_prime * Ez

        # Update position
        x += vx * dt
        y += vy * dt
        z += vz * dt

        # Update phase
        if !isnothing(wave)
            ω_t = wave_frequency(wave.modulation, wave, t_current)
            j = searchsortedfirst(wave.plasma.λ, λ_current)
            j = clamp(j, 1, length(wave.plasma.Ω_e))
            k_parallel = dispersion_relation_whistler_branch(ω_t, 0.0, wave.plasma.ω_pe[j], wave.plasma.Ω_e[j])
            v_parallel = vx*b̂x + vy*b̂y + vz*b̂z
            Ψ += k_parallel * v_parallel * dt - ω_t * dt
        end

        # Update stored values with current values
        t_current += dt

        # Update current s with the one closest to the grid
        if !isnothing(particle.field_line)
            dists = sqrt.((particle.field_line.xs_grid .- x).^2 .+
                          (particle.field_line.ys_grid .- y).^2 .+
                          (particle.field_line.zs_grid .- z).^2)
            s = particle.field_line.s_grid[argmin(dists)]
        end

        if store_trajectory
            r_array[i, :] .= (x, y, z)
            v_array[i, :] .= (vx, vy, vz)
            t_array[i] = t_current

            v_mag = sqrt(vx^2 + vy^2 + vz^2)
            v_par = vx*b̂x + vy*b̂y + vz*b̂z
            α_array[i] = acos(clamp(abs(v_par) / v_mag, 0.0, 1.0))
        end

        if x^2 + y^2 + z^2 ≤ r_ionosphere^2
            println("Precipitation!")
            if store_trajectory
                return (
                    position=r_array[1:i, :],
                    velocity=v_array[1:i, :],
                    time=t_array[1:i],
                    α=α_array[1:i],
                    footpoint=(x, y, z)
                )
            else
                return (time=t_current, footpoint=(x, y, z))
            end
        end
    end

    println("Bouncing...")
    if store_trajectory
        return (
            position=r_array,
            velocity=v_array,
            time=t_array,
            α=α_array,
            footpoint=nothing
        )
    else
        return (time=t_current, footpoint=nothing)
    end
end
