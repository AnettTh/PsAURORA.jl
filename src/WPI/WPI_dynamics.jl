using AURORA
using AURORA; mₑ, eV_in_J, c₀
using QuadGK
using LinearAlgebra
using Roots

#===

TODO: Several places in this file, the change in arc-length with latitude (ds/dλ) uses a
dipole approximation, also for the Tsyganenko model. This must either be fixed or jusitfied!

===#

"""
    wave_transit(
    ω,
    λ_resonance,
    particle::ParticleState,
    plasma::PlasmaState;
    field_dependent::Bool=true,
    )

Calculates the transit time of the whistler wave from the equator to the resonance latitude.

By either using a field-independent model (Miyoshi and Saito, 2012) or a field-dependent
model (Hsieh et al. 2022), where the latter is default. The function calculates the time it
takes for the wave specified by the frequency `ω`, the particle- and plasma state, traveling
from the source region (equator) to the resonance latitude `λ_resonance`. Note that the
result is given as an absolute value, meaning that there is no information about the
direction of the wave propagation.

# Arguments

- `ω`: Wave angular frequency of the whistler wave [rad/s].
- `λ_resonance`: The latitude of the resonance region (0.0 is equator) [rad].
- `particle`: Structure holding the particle state.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `field_dependent`: Decides which model to use, default is the field-dependent one.

# Returns

- Wave transit time for each frequency ω [s].
"""
function wave_transit(
    ω,
    λ_resonance,
    particle::ParticleState,
    plasma::PlasmaState;
    field_dependent::Bool=true,
    )

    # Initial position of the wave in the equatorial plane
    R0 = particle.L * RE

    # Define function to integrate
    function f(λ)
        # If field-dependent, choose based on position, else use equatorial Ω_e
        Ω_e = field_dependent ?
            Ωe_at_λ(λ, particle.L, particle.ϕ_eq, particle.magnetic_field) :
            plasma.Ω_e[1]

        ω_pe = field_dependent ?
            ωpe_at_λ(plasma.ne_model, particle.L, λ) :
            plasma.ω_pe[1]

        v_g = group_velocity_whistler_wave(ω, Ω_e, ω_pe)

        # Arc length element along a dipole field line
        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        return ds_dλ / v_g
    end

    # Integrate from equator to the resonance region and return result
    t_w, _ = quadgk(f, 0.0, λ_resonance)
    return abs.(t_w)
end


"""
    particle_transit(
    particle,
    λ_resonance;
    wave::Union{WaveState, Nothing}=nothing
    V_r=nothing,
    field_dependent::Bool=true,
    use_boris::Bool=false,
    )

Calculates the transit time of the particle from the resonance latitude to the ionosphere.

Either using a field-independent model (Miyoshi and Saito, 2012) or a field-dependent
model (Hsieh et al. 2022), the transit time, i.e. the time it takes a particle to travel
from the resonance region `λ_resonance` to the ionsophere, is calculated. The particle can
either travel freely, by setting `wave=nothing`, or can be subject to wave-particle
interactions, by passing also a `wave` to the function. To allow for the use of different
resonances `n`, the resonance velocity `V_r` can also be passed, allowing the function to
decide the hemisphere of the resonance latitude `λ_res` and direction of travel.

# Arguments

- `particle`: Structure holding the particle state.
- `λ_resonance`: The latitude of the resonance region (0.0 is equator) [rad].

# Keyword Arguments

- `wave`: Struct holding information about the wave (modulation and amplitude), default is
  `nothing`, i.e. no WPI considered.
- `V_r`: Resonance velocity of the particle [m/s].
- `field_dependent`: Decides which model to use, default is `true`.
- `use_boris`: Option to use the boris-mover to find the time-of-flight, default is `false`.

# Returns

- Particle transit time for each frequency ω [s].
"""
function particle_transit(
    particle,
    λ_resonance;
    wave::Union{WaveState, Nothing}=nothing,
    V_r=nothing,
    field_dependent::Bool=true,
    use_boris::Bool=false
)

    # If Boris-mover is to be used, shortcuts the function
    if use_boris
        # Find the initial position for the boris-mover, i.e. where the particle resonates
        result = boris_mover(
            particle;
            wave=wave,
            λ_start=λ_resonance,
            store_trajectory=false,
            stop_when_bouncing=true
        )

        # If the particle precipitates, return the TOF
        return isnothing(result.footpoint) ? nothing : result.time
    end

    # Decides which hemisphere the particle travels towards
    λ_iono = if isnothing(V_r)
        sign(λ_resonance) * particle.λ_ionosphere   # If non given, use λ at resonance
    else
        -sign(V_r) * particle.λ_ionosphere          # If given, base on velocity direction
    end

    # Function along a dipole fieldline to be integrated over
    function f(λ)
        if field_dependent
            # Velocity as a function of latitude
            B_λ = particle.magnetic_field(Spherical(particle.L, λ, particle.ϕ_eq))
            α_λ = pitch_angle_at_λ(particle.α_eq, particle.B_eq, norm(B_λ))

            # Find velocity, but removing mirrored particles
            isnothing(α_λ) && return 0.0
            vz = particle.v * cos(α_λ)
            iszero(vz) && return 0.0
        else
            # Velocity, but constant
            vz = abs(particle.v * cos(particle.α_lc))
        end

        # Position in the equatorial plane
        R0 = particle.L * RE

        # Arc length element along a dipole field line
        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        result = ds_dλ / vz

        # Removing NaN-values to not break the integration
        isnan(result) && return 0.0
        return result
    end

    # Integrate from the resonance region to the ionosphere and return
    t_e, _ = quadgk(f, λ_resonance, λ_iono)
    return abs.(t_e)
end


"""
    resonance_latitude(
    ω::Float64,
    particle::ParticleState,
    plasma::PlasmaState;
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false,
    wave_direction=-1
)

    resonance_latitude(
    ω_grid::AbstractVector,
    particle::ParticleState,
    plasma::PlasmaState;
    kwargs...
)

Calculate the resonance latitude for the defined particle with a whistler mode wave.

Depending on the frequency of the wave, the state of the particle and the state of the
plasma, the resonance latitude in a magnetic field is found. `n` denotes the harmonic, `1`
being the default value. Also returns the resonant velocity, used to decide where the
particles resonate from. For example, for n=1, the wave resonate with particles moving in
the opposite direction of the wave vector, meaning that particles precipitating in the
northern hemisphere comes from latitudes in the southern hemisphere.

# Arguments

- `ω` or `ω_grid`: Either the frequency or a grid of frequencoes of the wave [rad/s].
- `particle`: Structure holding the state of the particle.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `θ`: WNA, default is `0.0` (ducted).
- `n`: Resonance number, default is `1` (cyclotron resonance).
- `fast`: Option to use pre-computed values to also decide field strength. Avoids evaluating
  the Tsyganenko model (i.e. is faster), but result is less smooth. Default is `false`.
- `wave_direction`: Decides which direction the wave launces, default being `-1` which is
  anti-parallel with the magnetic field.

# Returns

- The latitude of resonance for the given frequency, particle and plasma [rad].
- The resonance velocity, used to decide which hemsiphere the resonating particles are
  located in [m/s].
"""
function resonance_latitude(
    ω::Float64,
    particle::ParticleState,
    plasma::PlasmaState;
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false,
    wave_direction=-1
)
    if abs(wave_direction) != 1
        throw(ArgumentError("Wave direction needs to be `-1` for southward or `1` for nothward."))
    end

    # Define the resonance condition (equation that should be zero for resonance to occure)
    function resonance_condition(λ; n=n, θ=θ, fast=fast)

        if fast
            # Choose latitude that matches closest to the precomputed grid
            i = searchsortedfirst(plasma.λ, λ)
            i = clamp(i, 1, length(plasma.Ω_e))
            Ω_e = plasma.Ω_e[i]
            ω_pe = plasma.ω_pe[i]

            # Compute field strength from known values
            B_mag = Ω_e * mₑ / abs(qₑ)
        else
            # Evaluate magnetic field at current position
            B = particle.magnetic_field(Spherical(particle.L, λ, particle.ϕ_eq))
            B_mag = norm(B)

            # Plasma-parameters based on current latitude
            Ω_e = abs(qₑ) * B_mag / mₑ
            ω_pe = ωpe_at_λ(plasma.ne_model, particle.L, λ)
        end

        # Find parallel component of wavenumber
        k = dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)
        k_parallel = wave_direction * k * cos(θ)

        # Find parallel velocity, keeping the sign of μ and returning NaN if mirrored
        α = pitch_angle_at_λ(particle.α_eq, particle.B_eq, B_mag)

        isnothing(α) && return NaN
        v_parallel = particle.v * cos(α) * (-sign(particle.μ))

        # Return resonance condition for given λ (from for example eq.4 in Hsieh 2022)
        return ω - k_parallel * v_parallel - (n*Ω_e / particle.γ)
    end

    # Find λ that causes the resonance condition to change sign, return NaN if there is none
    λ_grid = plasma.λ
    f_possible = resonance_condition.(λ_grid, n=n)

    idx = nothing

    # Look for index of change in sign
    for i in 1:length(f_possible)-1
        f_i = f_possible[i]
        f_i1 = f_possible[i+1]
        if !isnan(f_i) && !isnan(f_i1) && f_i * f_i1 < 0
            idx = i
            break
        end
    end

    # Return if none is found
    isnothing(idx) && return NaN, NaN

    # Figure out the latitude where the funciton changed sign, i.e. index `idx`
    λ_resonance = find_zero(resonance_condition, (λ_grid[idx], λ_grid[idx+1]))

    # Find also the resonance velocity at this position to decide the hemisphere
    i_res = searchsortedfirst(plasma.λ , λ_resonance)
    i_res = clamp(i_res, 1, length(plasma.Ω_e))
    Ω_e_res = plasma.Ω_e[i_res]
    ω_pe_res = plasma.ω_pe[i_res]
    k_res = dispersion_relation_whistler_branch(ω, θ, ω_pe_res, Ω_e_res)
    k_par_res = wave_direction * k_res * cos(θ)

    # Calculate the resonance velocity (Saito et al. 2012 #NOTE: Verify this!)
    V_r = (ω - (n * Ω_e_res / particle.γ)) / k_par_res

    return λ_resonance, V_r
end

function resonance_latitude(
    ω_grid::AbstractVector,
    particle::ParticleState,
    plasma::PlasmaState;
    kwargs...
)
    # Loop over the function if there is several frequencies
    λ_out = zeros(length(ω_grid))
    Vr_out = zeros(length(ω_grid))

    for (i, ω) in enumerate(ω_grid)
        result = resonance_latitude(ω, particle, plasma; kwargs...)

        if isnan(result[1])
            λ_out[i]  = NaN
            Vr_out[i] = NaN
        else
            λ_out[i]  = result[1]
            Vr_out[i] = result[2]
        end
    end

    return λ_out, Vr_out
end


"""
    WPI_TOF(
    ω_grid::AbstractVector,
    particle::ParticleState,
    plasma::PlasmaState;
    wave::Union{WaveState, Nothing}=nothing,
    use_boris::Bool=false,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false
)

    WPI_TOF(ω::Real, particle::ParticleState, plasma::PlasmaState; kwargs...)

Calculates the time-of-flight (TOF) of a particle resonating with a whistler mode chorus
wave, then precipitating into the ionosphere.

The function uses a defined frequency-grid for the wave, as well as a struct holding the
particle state and the plasma state, to calculate the TOF. If a `wave` is passed to the
function, also wave-particle-interactions (WPI) are considered. The TOF can be calculated
either by integration (using field-dependent or -independent models), or by numerical
simulations (boris-mover), the latter by setting `use_boris = true`. Also the resonance `n`
can be decided, default is `1` (cyclotron resonance), and the wave-normal-angle `θ`, default
is `0.0`. The function can be run slightly faster by using pre-computed values for the
electron plasma frequency `ω_pe` and the electron cyclotron frequency `Ω_e`, by setting
`fast=true` (default is `false`), but this gives less smooth results.

# Arguments

- `ω_grid`: The frequency grid of the wave [rad/s].
- `particle`: Structure holding the state of the particle.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `wave`: Struct holding information about the wave to consider (if there is one), default
  is `nothing`.
- `use_boris`: Decides whether to integrate or use the boris-mover, default is `true`.
- `field_dependent`: Decides which model to use, default is `true`.
- `θ`: Wave-normal angle of the wave, default is `1.0` (field-aligned).
- `n`: Harmonic number, default is `1` (cyclotron resonance).
- `fast`: Option to use precomputed values for the PlasmaState, default is `false`.

# Returns

- The time of flight for the `particle` in the current state of the `plasma` under the
  influence of the defined `wave`, using the defined method [s]. If the particle don't
  precipitate, the TOF returned is `0.0`.
"""
function WPI_TOF(
    ω_grid::AbstractVector,
    particle::ParticleState,
    plasma::PlasmaState;
    wave::Union{WaveState, Nothing}=nothing,
    use_boris::Bool=false,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false
)

    # Calculate wave-launch time based on given model
    t_l = wave_launch_time(wave, ω_grid)

    tof = zeros(length(ω_grid))

    # Loop over all frequencies as they have different resonance regions
    for (i, ω) in enumerate(ω_grid)
        λ_res, V_r = resonance_latitude(ω, particle, plasma; θ=θ, n=n, fast=fast)

        # For conbinations that don't resonate
        if isnothing(λ_res) || isnan(λ_res)
            tof[i] = NaN
            continue
        end

        # Calculate the components for the given ω and sum them
        t_w = wave_transit(
            ω,
            λ_res,
            particle,
            plasma;
            field_dependent=field_dependent
        )

        t_e = particle_transit(
            particle,
            λ_res;
            wave=wave,
            V_r=V_r,
            use_boris=use_boris,
            field_dependent=field_dependent,
        )

        # Total time = wave-to-λ + particle-to-ionosphere + delay of current ω based on wave
        if !isnothing(t_e)
            tof[i] = t_w + t_e + t_l[i]
        else
            tof[i] = 0.0  # TODO: Is there a better way to do this maybe?
        end
    end
    return tof
end

# Also allows for single frequencies
function WPI_TOF(ω::Real, particle::ParticleState, plasma::PlasmaState; kwargs...)
    return WPI_TOF([ω], particle, plasma; kwargs...)[1]
end
