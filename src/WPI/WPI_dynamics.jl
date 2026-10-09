using AURORA
using AURORA; mₑ, eV_in_J, c₀
using QuadGK
using LinearAlgebra
using Roots

#===
TODO: Several places in this file, the change in arc-length with latitude (ds/dλ) uses a
dipole approximation, also for the Tsyganenko model. This must either be fixed or jusitfied!

NOTE: When computing the resonance latitude, need to make sure that the wave propagation can
be controlled. Right now, it is tested for waves propagating south, such that the resonance
of particles precipitating in the northern hemisphere will resonate in the southern
hemisphere, but want to let also the opposite be possible.
===#
"""
    wave_transit(
    ω,
    λ_resonance,
    particle::ParticleState,
    plasma::PlasmaState;
    field_dependent::Bool=true)

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
    field_dependent::Bool=true
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


# TODO: Make sure this traces from the southern resonance region to the norhtern hemisphere!
"""
    particle_transit(particle, λ_resonance; V_r=nothing z_ionosphere::Float64=600e3)

Calculates the transit time of the particle from the resonance latitude to the ionosphere.

By either using a field-independent model (Miyoshi and Saito, 2012) or a field-dependent
model (Hsieh et al. 2022), where the latter is default. The function calculates the time it
takes for the particle in the defined state to travel from the resonance region
`λ_resonance` to the defined ionospheric latitude, `λ_ionosphere`, as defined in the
particle state (600 km altitude).

# Arguments

- `particle`: Structure holding the particle state.
- `λ_resonance`: The latitude of the resonance region (0.0 is equator) [rad].

# Keyword Arguments

- `field_dependent`: Decides which model to use, default is `true`.
- `use_boris`: Option to use the boris-mover to find the time-of-flight, default is `false`.

# Returns

- Particle transit time for each frequency ω [s].
"""
function particle_transit(
    particle,
    λ_resonance;
    V_r=nothing,
    field_dependent::Bool=true,
    use_boris::Bool=false
)

    # If Boris-mover is to be used, shortcuts the function
    if use_boris
        # Find the initial position for the boris-mover, i.e. where the particle resonates
        result = boris_mover(particle; λ_start=λ_resonance, store_trajectory=false)
        return result.time
    end

    # Position  in the equatorial plane
    R0 = particle.L * RE

    # Decides which hemisphere the particle travels towards
    λ_iono = if isnothing(V_r)
        sign(λ_resonance) * particle.λ_ionosphere   # If non given, use λ at resonance
    else
        -sign(V_r) * particle.λ_ionosphere          # If given, base on velocity direction
    end

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

        # Arc length element along a dipole field line
        ds_dλ = R0 * sqrt(1 + 3sin(λ)^2) * cos(λ)
        result = ds_dλ / vz

        # Removing NaN-values to not break the integration
        isnan(result) && return 0.0
        return result
    end

    # Integrate from the resonance region to the ionosphere and return
    t_e, _ = quadgk(f, λ_resonance, λ_iono)  # NOTE: takes opposite hemisphere as default
    return abs.(t_e)
end


# TODO: Use AbstractWave instead of this
"""
    wave_chirp(ω::AbstractVector; t::Real=0.2)

Compute the chorus wave launch time profile t₀(ω) at the magnetic equator.

Returns the time at which a wave of angular frequency ω is launched from the equatorial
source region, assuming a linear frequency chirp (rising tone), as modeled in Chen 2020. The
default parameter of rise-time in 0.2 s correspond to Chen 2020, and the start- and end
frequencies is the start- and end point of the frequency range.

# Arguments

- `ω`: Angular frequency at which to evaluate t₀ [rad/s].

# Keyword Arguments

- `t`: Duration of the chorus element [s], default is 0.2 s, similar to what is used in Chen
  2021.

# Returns

- A function for release time of each ω.
"""
function wave_chirp(ω::AbstractVector; t::Real=0.2)
    chirp_rate = (ω[1] - ω[end]) / t
    return (ω .- ω[1]) ./ chirp_rate
end


"""
    resonance_latitude(
    ω::Float64,
    particle::ParticleState,
    plasma::PlasmaState;
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false
)

    resonance_latitude(
    ω_grid::AbstractVector,
    particle,
    plasma;
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false
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

- `θ`: WNA, default is nothing (ducted).
- `n`: Resonance number, default is `1` (cyclotron resonance).
- `fast`: Option to use pre-computed values to also decide field strength. Avoids evaluating
  the Tsyganenko model (i.e. is faster), but result is less smooth. Default is `false`.

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
    wave_direction=-1       # NOTE: South???
)
    if abs(wave_direction) != 1
        throw(ArgumentError("Wave direction needs to be `-1` for southward or `1` for nothward."))
    end

    # Define the resonance condition (equation that should equal zero)
    function resonance_condition(λ; n=n, θ=θ, fast=fast)   # TODO: Add full spherical coordinates as an argument?

        if fast
            # Choose latitude that matches closest to the precomputed grid
            i   = searchsortedfirst(plasma.λ, λ)
            i   = clamp(i, 1, length(plasma.Ω_e))
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

    # Find which λ causes the resonance condition-function to change sign
    λ_grid = plasma.λ
    f_possible = resonance_condition.(λ_grid, n=n)

    # Find which λ causes the function to change sign, return NaN if there is none
    idx = nothing
    for i in 1:length(f_possible)-1
        if !isnan(f_possible[i]) && !isnan(f_possible[i+1]) && f_possible[i] * f_possible[i+1] < 0
            idx = i
            break
        end
    end
    isnothing(idx) && return NaN, NaN

    # Figure out the latitude where the funciton changed sign, i.e. index `idx`
    λ_resonance = find_zero(resonance_condition, (λ_grid[idx], λ_grid[idx+1]))

    # Find also the resonance velocity at this position to decide the hemisphere of the
    # resonant particles
    i_res = searchsortedfirst(plasma.λ , λ_resonance)
    i_res = clamp(i_res, 1, length(plasma.Ω_e))
    Ω_e_res = plasma.Ω_e[i_res]
    ω_pe_res = plasma.ω_pe[i_res]
    k_res = dispersion_relation_whistler_branch(ω, θ, ω_pe_res, Ω_e_res)
    k_par_res = wave_direction * k_res * cos(θ)
    V_r = (ω - (n * Ω_e_res / particle.γ)) / k_par_res

    return λ_resonance, V_r
end

function resonance_latitude(
    ω_grid::AbstractVector,
    particle,
    plasma;
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false
)
    # Loop over the function if there is several frequencies
    λ_out = zeros(length(ω_grid))
    Vr_out = zeros(length(ω_grid))

    for (i, ω) in enumerate(ω_grid)
        result = resonance_latitude(ω, particle, plasma; θ=θ, n=n, fast=fast)

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
    use_boris::Bool=false,
    wave_launch_time=wave_chirp,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1
)

    WPI_TOF(ω::Real, particle, plasma; kwargs...)

Calculates the time-of-flight of a particle resonating with a whistler mode chorus wave,
then precipitating into the ionosphere.

The funciton uses a defined frequency-grid for the wave, particle state and plasma state to
calculate the resonance latitude. The wave launch time can be given as a function in the
keyword arguments, else, it is assumed instant at all frequencies.

# Arguments

- `ω_grid`: The frequency grid of the wave [rad/s].
- `particle`: Structure holding the state of the particle.
- `plasma`: Structure holding the state of the plasma on a latitude grid.

# Keyword Arguments

- `use_boris`: Decides whether to integrate or use the boris-mover, default is `true`.
- `wave_launch_time`: The time at which the frequencies in the frequency grid is launched
  from the source region (equator). The default is a simple chirp, but can be any function
  that takes frequencies `ω`.
- `field_dependent`: Decides which model to use, default is `true`.
- `θ`: Wave-normal angle of the wave, default is `1.0` (field-aligned).
- `n`: Harmonic number, default is `1` (cyclotron resonance).

# Returns

- The time of flight for the `particle` in the current state of the `plasma` under the
  influence of the defined wave `wave_launch_time` [s].
"""
function WPI_TOF(
    ω_grid::AbstractVector,
    particle::ParticleState,
    plasma::PlasmaState;
    use_boris::Bool=false,
    wave_launch_time=wave_chirp,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1,
    fast=false
)

    # Calculate wave-launch time based on given model
    # NOTE: There should be kwargs here, but will have to make this into a struct anyways
    t_l = wave_launch_time(ω_grid)

    tof = zeros(length(ω_grid))

    # Precompute magnetic-field before frequency loop
    if use_boris
        if isnothing(particle.B_interpolated)
            # TODO: remove this when non-precomputed B-field is removed
            #throw(ArgumentError("Trying to remove non-precomputed magnetic field..."))
            B_interp = get_magnetic_field(particle.magnetic_field, particle.L, particle.ϕ_eq)
        else
            B_interp = particle.B_interpolated
        end
    else
        B_interp = nothing
    end

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
            V_r=V_r,
            use_boris=use_boris,
            field_dependent=field_dependent,
        )

        # Total time = wave-to-λ + particle-to-ionosphere + delay of current ω based on wave
        tof[i] = t_w + t_e + t_l[i]

    end
    return tof
end

# Also allows for single frequencies
function WPI_TOF(ω::Real, particle, plasma; kwargs...)
    return WPI_TOF([ω], particle, plasma; kwargs...)[1]
end
