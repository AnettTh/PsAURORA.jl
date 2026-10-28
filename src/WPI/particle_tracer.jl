using AURORA
using ProgressMeter

"""
    boris_mover(
    particle::ParticleState;
    wave::Union{WaveState, Nothing}=nothing,
    λ_start=nothing,
    φ0::Float64=0.0,
    n_T::Int=100_000,
    resolution::Int=10,
    store_trajectory::Bool=false
    count_bounces::Bool=false,
    stop_when_bouncing::Bool=false,
)

Simulate the trajectory of a charged particle using the Boris mover.

Uses the information about the particle given by `ParticleState` to numerically trace its
gyromotion in a magnetic field. The function also allows for wave-particle interaction (WPI)
by adding perturbations to the magnetic- and electric fields. Depending on the purpose of
the simulation, it can be run while storing trajectory, counting bounces or stopping when a
bounce is detected.

# Arguments

- `particle`: Struct holding the state of the particle.

# Keyword Arguments

- `wave`: Struct holding information about the wave, default is `nothing`, which does not
  consider WPI.
- `λ_start`: Initial latitude to start the particle from, allows simulation to start at the
  position of resonance between the particle and wave [rad], default is `nothing`, which
  uses the particles initial position to calculate the initial latitude.
- `φ0`: Initial phase of the velocity vector for the particle, default is `0.0` [rad].
- `n_T`: Number of (initial) gyroperiods to simulate for, default is `100_000`.
- `resolution`: Steps per gyroperiod, default is `10`.
- `store_trajectory`: Option to store the trajectory, default is `false`.
- `count_bounces`: Option to track how many bounces are done before timing out (or
  precipitating), default is `false`.
- `stop_when_bouncing`: Option to stop the simulation if the particle bounces, default is
  `false`.

# Returns

- If the trajectroy is stored, a named tuple containing `position`, `velocity`, `time`, `α`,
  `footpoint`, if not, only the time-of-flight `time` and the position of impact with the
  ionsophere, `footpoint` is returned. Also returns number of bounces conducted, if
  counting, else this is just an empty variable.

# Throws

- `ArgumentError`: If the maximum number of gyrophases or the number of samples per period
  is set to be a negative integer.
"""
function boris_mover(
    particle::ParticleState;
    wave::Union{WaveState, Nothing}=nothing,
    λ_start::Union{Float64, Nothing}=nothing,
    φ0::Float64=0.0,
    n_T::Int=100_000,
    resolution::Int=10,
    store_trajectory::Bool=false,
    count_bounces::Bool=false,
    stop_when_bouncing::Bool=false,
)

    #===

    Check that the input is valid before starting the simulation.

    ===#
    n_T > 0 || throw(ArgumentError("n_T must be positive"))
    resolution > 0 || throw(ArgumentError("resolution must be positive"))

    #===

    Initialize the particle.

    ===#

    # Define where the particle has hit the ionosphere
    r_ionosphere = RE + z_ionosphere

    # Find total number of steps reqired
    steps = Int(n_T * resolution)

    if store_trajectory
        @info "\nStoring trajectory"
        r_array = zeros(steps + 1, 3)
        v_array = zeros(steps + 1, 3)
        t_array = zeros(steps + 1)
        α_array = zeros(steps + 1)
    end

    # Find position and velocity
    if isnothing(λ_start)
        @info "\nUsing equatorial position as initial position"
        # Use default positions, usually in the equatorial plane
        x, y, z = particle.r0
        vx, vy, vz = particle.v0
    else
        @info "\nUsing latitude λ=$(round(rad2deg(λ_start), digits=1))° as initial position"
        # Get position of resonance
        # TODO: Using ϕ_eq here is just temporary, figure out the accurate position?
        r_start = Cartesian(Spherical(particle.L, λ_start, particle.ϕ_eq))
        x, y, z = r_start.x, r_start.y, r_start.z

        # Find v0 at resonance latitude: Evaluate the exact magnetic field once
        B_start = particle.magnetic_field(Spherical(particle.L, λ_start, particle.ϕ_eq))
        b̂, e1, e2, B_start_mag = magnetic_basis(B_start)

        # Use initial pitch angle to decompose velocity
        α_start = pitch_angle_at_λ(particle.α_eq, particle.B_eq, B_start_mag)
        v_par = -particle.v * cos(α_start) * sign(particle.μ)
        v_perp = particle.v * sin(α_start)

        vx, vy, vz = v_par .* b̂ .+ v_perp * cos(φ0) .* e1 .+ v_perp * sin(φ0) .* e2
    end

    # Variable for tracking the time
    t_current = 0.0

    # Variable to detect and count bounces
    if count_bounces || stop_when_bouncing
        n_bounces = 0
        v_par_previous = vx*b̂[1] + vy*b̂[2] + vz*b̂[3]
    end

    # Store initial values
    if store_trajectory
        r_array[1, :] .= (x, y, z)
        v_array[1, :] .= (vx, vy, vz)
        t_array[1] = t_current
        α_array[1] = acos(abs(particle.μ))
    end

    # Keep track of phase if wave is to be considered
    Ψ = isnothing(wave) ? 0.0 : wave.Ψ0     # TODO: change name to Ψ_current, for consistency?

    #===

    Update the particle and fill the arrays containing the trajectory

    ===#
    for i in 2:(steps + 1)

        # Background fields
        B = particle.magnetic_field(Cartesian(x, y, z))
        Bx, By, Bz = B[1], B[2], B[3]
        Ex, Ey, Ez = particle.E0

        # Add wave perturbation at current time-step, if a wave exists
        if !isnothing(wave)
            # TODO: This is a dipole approximation, is there a better way?
            r_mag = sqrt(x^2 + y^2 + z^2)
            λ_current = asin(z / r_mag)

            # TODO: Placeholder, this only matters if the amplitude of the wave depends on the position
            s_current = 0.0

            # Find perturbation and update the current field
            δB, δE = wave_fields(
                wave,
                t_current,
                s_current,
                Ψ,
                SVector(Bx, By, Bz),
                λ_current
            )

            Bx += δB[1]
            By += δB[2]
            Bz += δB[3]

            Ex += δE[1]
            Ey += δE[2]
            Ez += δE[3]
        end

        # Find current magnitude and direction of the magnetic field
        B_mag = sqrt(Bx^2 + By^2 + Bz^2)
        b̂x, b̂y, b̂z = Bx/B_mag, By/B_mag, Bz/B_mag

        # Define cyclotron frequency and period for current timestep
        Ω_e = (B_mag * abs(qₑ)) / mₑ
        T = 2π / Ω_e

        # Define adaptive timestep
        dt = T / resolution

        # Half electric acceleration
        q_prime = dt * (-qₑ) / 2mₑ      # NOTE: qₑ here is signed, should this be done in the constant-file instead? Check consistency!

        #===

        This is the core Boris-scheme!

        ===#

        # Half electric acceleration
        vx_minus = vx + q_prime * Ex
        vy_minus = vy + q_prime * Ey
        vz_minus = vz + q_prime * Ez

        # Magnetic rotation vector
        hx = q_prime * Bx
        hy = q_prime * By
        hz = q_prime * Bz
        h_mag2 = hx^2 + hy^2 + hz^2

        # Some intermediate calculations...
        sx = 2hx / (1 + h_mag2)
        sy = 2hy / (1 + h_mag2)
        sz = 2hz / (1 + h_mag2)

        ux_prime = vx_minus + (vy_minus * hz - vz_minus * hy)
        uy_prime = vy_minus + (vz_minus * hx - vx_minus * hz)
        uz_prime = vz_minus + (vx_minus * hy - vy_minus * hx)

        vx_plus = vx_minus + (uy_prime * sz - uz_prime * sy)
        vy_plus = vy_minus + (uz_prime * sx - ux_prime * sz)
        vz_plus = vz_minus + (ux_prime * sy - uy_prime * sx)

        # Update the velocity
        vx = vx_plus + q_prime * Ex
        vy = vy_plus + q_prime * Ey
        vz = vz_plus + q_prime * Ez

        # Check if the particle has bounced in this time-step
        if stop_when_bouncing || count_bounces
            v_par_current = vx*b̂x + vy*b̂y + vz*b̂z

            # If it bounces, either count it or return it, depending on the config
            if sign(v_par_previous) != sign(v_par_current)

                if count_bounces
                    n_bounces += 1
                end

                if stop_when_bouncing
                    @info "\nBounced (and told to stop at that), no footpoint detected"
                    if store_trajectory
                        return (
                            position=r_array[1:i, :],
                            velocity=v_array[1:i, :],
                            time=t_array[1:i],
                            α=α_array[1:i],
                            footpoint=nothing,
                            n_bounces=count_bounces ? n_bounces : nothing,
                        )
                    else
                        return (
                            time=t_current,
                            footpoint=nothing,
                            n_bounces=count_bounces ? n_bounces : nothing,
                        )
                    end
                end
            end
            # If not, update current parallel velocity to previous parallel velocity
            v_par_previous = v_par_current
        end

        # Update the position
        x += vx * dt
        y += vy * dt
        z += vz * dt

        # Update the phase, if wave is to be considered
        if !isnothing(wave)
            # Find current frequency (based on modulation)
            ω_t = wave_frequency(wave, t_current)

            # Find current wave vector based on current plasma state and wave frequency
            j = searchsortedfirst(wave.plasma.λ, λ_current)
            j = clamp(j, 1, length(wave.plasma.Ω_e))
            θ = 0.0     # NOTE: Does WNA need to be 0 for all possible models?
            k_parallel = dispersion_relation_whistler_branch(
                ω_t,
                θ,
                wave.plasma.ω_pe[j],
                wave.plasma.Ω_e[j]
            )

            # The current velocity parallel to the magnetic field b̂
            v_parallel = vx*b̂x + vy*b̂y + vz*b̂z

            # Update the wave-phase (based on Saito et al. 2012)
            Ψ += k_parallel * v_parallel * dt - ω_t * dt
        end

        # Update stored values with current values
        t_current += dt

        if store_trajectory
            r_array[i, :] .= (x, y, z)
            v_array[i, :] .= (vx, vy, vz)
            t_array[i] = t_current

            v_mag = sqrt(vx^2 + vy^2 + vz^2)
            v_par = vx*b̂x + vy*b̂y + vz*b̂z
            α_array[i] = acos(clamp(abs(v_par) / v_mag, 0.0, 1.0))
        end


        #===

        Return if the particle precipitates

        ===#
        if x^2 + y^2 + z^2 ≤ r_ionosphere^2
            @info "\nPrecipitation!"
            if store_trajectory
                return (
                    position=r_array[1:i, :],
                    velocity=v_array[1:i, :],
                    time=t_array[1:i],
                    α=α_array[1:i],
                    footpoint=(x, y, z),
                    n_bounces=count_bounces ? n_bounces : nothing,
                )
            else
                return (
                    time=t_current,
                    footpoint=(x, y, z),
                    n_bounces=count_bounces ? n_bounces : nothing,
                )
            end
        end
    end

    #===

    Return result if the simulation times out (n_T is exceeded).

    ===#

    @info "\nBouncing..."
    if store_trajectory
        return (
            position=r_array,
            velocity=v_array,
            time=t_array,
            α=α_array,
            footpoint=nothing,
            n_bounces=count_bounces ? n_bounces : nothing,
        )
    else
        return (
            time=t_current,
            footpoint=nothing,
            n_bounces=count_bounces ? n_bounces : nothing,
        )
    end
end
