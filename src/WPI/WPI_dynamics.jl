using AURORA
using AURORA; mₑ, eV_in_J, c₀
using QuadGK
using LinearAlgebra
using Roots

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
    fieldline::AbstractFieldLine,
    plasma::PlasmaState;
    field_dependent::Bool=true,
    )

    # Index of the resonance and equator on the line
    i_res = clamp(searchsortedfirst(fieldline.λs, λ_resonance), 1, length(fieldline.λs))
    i_eq = fieldline.i_eq

    # Get segment of the line where the wave travels
    i_wave = min(i_eq, i_res):max(i_eq, i_res)
    s_wave = fieldline.s[i_wave]

    if field_dependent
        Ω_e = plasma.Ω_e[i_wave]
        ω_pe = plasma.ω_pe[i_wave]
    else
        # Use the values at the magnetic equator
        Ω_e  = fill(plasma.Ω_e[i_eq], length(i_wave))
        ω_pe = fill(plasma.ω_pe[i_eq], length(i_wave))
    end

    # The travel-time is found using the trapezoid method
    inverse_vg = @. 1 / group_velocity_whistler_wave(ω, Ω_e, ω_pe)
    t_w = sum(0.5 .* (inverse_vg[1:end-1] .+ inverse_vg[2:end]) .* abs.(diff(s_wave)))
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
    particle::ParticleState,
    fieldline::AbstractFieldLine,
    λ_resonance;
    wave::Union{WaveState, Nothing}=nothing,
    field_dependent::Bool=true,
    use_boris::Bool=false,
    get_footpoint::Bool=false,
)

    # Index on the line where the particle and wave meets
    i_res = clamp(searchsortedfirst(fieldline.λs, λ_resonance), 1, length(fieldline.λs))

    # If Boris-mover is to be used, shortcuts the function
    if use_boris
        r_resonance = Cartesian(fieldline.xs[i_res], fieldline.ys[i_res], fieldline.zs[i_res])
        # Find the initial position for the boris-mover, i.e. where the particle resonates
        result = boris_mover(
            particle;
            wave=wave,
            store_trajectory=false,
            stop_when_bouncing=true,
            r_resonance=r_resonance,
        )

        # If the particle precipitates, return the TOF
        # TODO: Figure out how to solve this! One for bounce, one for "did not reach"?
        if !get_footpoint
            return isnothing(result.footpoint) ? nothing : result.time
        else
            return isnothing(result.footpoint) ?
            (time=nothing, footpoint=nothing) :
            (time=result.time, footpoint=result.footpoint)
        end
    end
    # Check if the particle moves in the correct direction
    toward_end = fieldline.B_along_index == (particle.μ0 < 0)
    i_end = toward_end ? length(fieldline.λs) : 1
    position_fp = toward_end ? fieldline.position_fp_north : fieldline.position_fp_south

    # Segment the particle travels, i.e. from the resonance to the ionosphere
    i_particle = i_res ≤ i_end ? (i_res:i_end) : (i_res:-1:i_end)
    s_particle = fieldline.s[i_particle]

    # Find the parallel speed at  each stored point
    if field_dependent
        # TODO: Check what this returns for a mirrored particle and match it up with the behavior
        α = [pitch_angle(particle.α_eq, particle.B_eq, B) for B in fieldline.Bs[i_particle]]
        # Check if the particle have mirrored, then return nothing
        any(isnothing, α) && return get_footpoint ? (time=nothing, footpoint=nothing) : nothing
        v_par = [particle.v0_mag * cos(a) for a in α]
    else
        v_par = fill(abs(particle.v0_mag * cos(particle.α_lc)), length(i_particle))
    end

    # Using the trapezoid method to integrate
    inverse_v = 1 ./ v_par
    t_e = sum(0.5 .* (inverse_v[1:end-1] .+ inverse_v[2:end]) .* abs.(diff(s_particle)))
    return get_footpoint ?
        (time=t_e, footpoint=(position_fp.x, position_fp.y, position_fp.z)) :
        t_e
end

# TODO: Update
"""
    point_of_resonance(
    ω::Float64,
    particle::ParticleState,
    plasma::PlasmaState;
    θ::Float64=0.0,
    n::Int=1,
    fast::Bool=false,
    wave_direction=-1
)

    point_of_resonance(
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
function point_of_resonance(
    ω::Float64,
    particle::ParticleState,
    fieldline::AbstractFieldLine,
    plasma::PlasmaState;
    θ::Float64=0.0,
    n::Int=1,
    wave_direction::Int=-1,
)
    # Add some throws
    abs(wave_direction) == 1 || throw(
        ArgumentError("Wave direction needs to be `-1` for southward or `1` for northward.")
    )
    length(plasma.λ) == length(fieldline.λs) || throw(
        ArgumentError("The plasma was not built from this fieldline")
    )

    # Decide on the direction of the wave and the particle
    wave_along_b = wave_direction * (fieldline.B_along_index ? 1 : -1)
    v_par_sign   = -sign(particle.μ0)

    # Calculate the resonance condition at point `i`, where zero gives resonance
    function resonance_condition(i)
        # Current frequencies
        Ω_e = plasma.Ω_e[i]
        ω_pe = plasma.ω_pe[i]

        # Check if the particle have mirrored already
        α = pitch_angle(particle.α_eq, particle.B_eq, fieldline.Bs[i])
        isnothing(α) && return NaN

        # Use the dispersion relation to find k, if any
        k = dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)

        # Pick out only parallel components
        k_parallel = wave_along_b * k * cos(θ)
        v_parallel = particle.v0_mag * cos(α) * v_par_sign

        # Return the equation that is to be zero for resonance to happen
        return ω - (k_parallel * v_parallel) - (n * Ω_e / particle.γ)
    end

    # Look at the resonance condition along the entire fieldline
    f = [resonance_condition(i) for i in eachindex(fieldline.λs)]

    # Find the indices where `f` changes sign and is not-empty
    is_crossing(j) = !isnan(f[j]) && !isnan(f[j+1]) && f[j] * f[j+1] < 0
    i = findfirst(is_crossing, 1:length(f)-1)
    isnothing(i) && return nothing

    # Use linear interpolation to find the root between point i and i+1
    weight = f[i] / (f[i] - f[i+1])
    intrp(value) = value[i] + weight * (value[i+1] - value[i])

    # Get the values to return
    Ω_e_resonance = intrp(plasma.Ω_e)
    ω_pe_resonance = intrp(plasma.ω_pe)
    λ_resonance = intrp(fieldline.λs)

    r_resonance = Cartesian(intrp(fieldline.xs), intrp(fieldline.ys), intrp(fieldline.zs))

    k_resonance = dispersion_relation_whistler_branch(ω, θ, ω_pe_resonance, Ω_e_resonance)
    v_resonance = (ω - (n*Ω_e_resonance/particle.γ)) / (wave_along_b * k_resonance * cos(θ))

    return (λ=λ_resonance, r=r_resonance, v=v_resonance, i=i)
end

# TODO: This can be done faster, as Ω_e, ω_pe and α don't change with ω
function point_of_resonance(
    ω_grid::AbstractVector,
    particle::ParticleState,
    fieldline::AbstractFieldLine,
    plasma::PlasmaState;
    kwargs...
)
    return [point_of_resonance(ω, particle, fieldline, plasma; kwargs...) for ω in ω_grid]
end

# TODO: Update docstring
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
    fieldline::AbstractFieldLine,
    plasma::PlasmaState;
    wave::Union{WaveState, Nothing}=nothing,
    use_boris::Bool=false,
    field_dependent::Bool=true,
    θ::Float64=0.0,
    n::Int=1,
    get_footpoint::Bool=false,
)

    # Calculate wave-launch time based on given model
    t_l = wave_launch_time(wave, ω_grid)

    tof = zeros(length(ω_grid))
    footpoint = fill(NaN, 3, length(ω_grid))

    # Loop over all frequencies as they have different resonance regions
    for (i, ω) in enumerate(ω_grid)
        resonance = point_of_resonance(ω, particle, fieldline, plasma; θ=θ, n=n)

        # For conbinations that don't resonate
        if isnothing(resonance) #|| isnan(λ_res)
            tof[i] = NaN
            continue
        end

        # Calculate the components for the given ω and sum them
        t_w = wave_transit(
            ω,
            resonance.λ,
            fieldline,
            plasma;
            field_dependent=field_dependent,
        )

        p_t_result = particle_transit(
            particle,
            fieldline,
            resonance.λ;
            wave=wave,
            use_boris=use_boris,
            field_dependent=field_dependent,
            get_footpoint=get_footpoint,
        )

        if get_footpoint
            t_e = isnothing(p_t_result) ? nothing : p_t_result.time
            fp  = isnothing(p_t_result) ? nothing : p_t_result.footpoint
            isnothing(fp) || (footpoint[:, i] .= fp)
        else
            t_e = p_t_result
        end

        # Total time = wave-to-λ + particle-to-ionosphere + delay of current ω based on wave
        tof[i] = isnothing(t_e) ? 0.0 : t_w + t_e + t_l[i]
    end
    return get_footpoint ? (tof=tof, footpoint=footpoint) : tof
end

# Also allows for single frequencies
function WPI_TOF(
    ω::Real,
    particle::ParticleState,
    fieldline::AbstractFieldLine,
    plasma::PlasmaState;
    kwargs...
)
    return WPI_TOF([ω], particle, fieldline, plasma; kwargs...)[1]
end
