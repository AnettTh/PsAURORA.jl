using AURORA
using StaticArrays

# Models for defining waves intearcting with the particles.
# Implements AbstractWaveModulation with MonochromaticModulation, ChirpModulation and Demekhov and
# AbstractWaveAmplitude with ConstantAmplitude and GaussianAmplitude (TBW).

"""
    AbstractWaveModulation

Abstract type for wave modulation models.
"""
abstract type AbstractWaveModulation end

"""
    AbstractWaveAmplitude

Abstract type for wave amplitude models.
"""
abstract type AbstractWaveAmplitude end

"""
    MonochromaticModulation

Monochromatic wave modulation where the frequency stays constant.

# Fields

- `ω`: Wave frequency [rad/s].
"""
struct MonochromaticModulation <: AbstractWaveModulation
    ω::Float64
end

"""
    ChirpModulation

Chirp wave modulation with a linear increase/decrease in frequency over some time.

# Fields

- `ω0`: Frequency at the start [rad/s].
- `ω1`: Frequency at the end [rad/s].
- `t_duration`: Total duration of chirp [s].
"""
struct ChirpModulation <: AbstractWaveModulation  # TODO: Simplify how this is called?
    ω0         :: Float64
    ω1         :: Float64
    t_duration :: Float64
end

# Chirp wave as described by Demekhov 2016 (TBW)
struct DemekhovModulation <: AbstractWaveModulation end


"""
    ConstantAmplitude

Constant amplitude model, which stays constant for all times `t` and positions `s`.

# Fields

- `Bw`: Amplitude of the wave [T].
"""
struct ConstantAmplitude <: AbstractWaveAmplitude
    Bw::Float64
end

wave_amplitude(a::ConstantAmplitude, t, s) = a.Bw

# TODO: Add also Gaussian-like amplitude from Chen 2020. Then, make sure to look into `s_current` in boris-mover


"""
    WaveState{M<:AbstractWaveModulation, A<:AbstractWaveAmplitude}

Collecting modulation and amplitude to define a complete wave state.

# Fields

- `modulation`: Struct holding information about the modulation.
- `amplitude`: Struct holding information about the amplitude.
- `plasma`: Struct holding information about the plasma state, used for dispersion relation.
- `Ψ0`: Initial wave phase [rad].
- `wave_direction`: Decides the propagation direction, `1` is northward and `-1` southward.

Can be constructed like this, for example:
    WaveState(
    modulation::AbstractWaveModulation,
    amplitude:AbstractWaveAmplitude,
    plasma::PlasmaState;
    Ψ0=2π*rand(),
    wave_direction=-1
    )
"""
struct WaveState{M<:AbstractWaveModulation, A<:AbstractWaveAmplitude}
    modulation     :: M
    amplitude      :: A
    plasma         :: PlasmaState
    Ψ0             :: Float64
    wave_direction :: Int
    wave_along_b   :: Int
end


function WaveState(
    modulation,
    amplitude,
    plasma,
    fieldline::AbstractFieldLine;
    Ψ0=2π*rand(),
    wave_direction=-1
)
    abs(wave_direction) == 1 || throw(ArgumentError("wave_direction must be 1 or -1"))

    length(plasma.λ_grid) == length(fieldline.λs) || throw(
        ArgumentError("The plasma was not built from this fieldline")
    )

    wave_along_b = wave_direction * (fieldline.B_along_index ? 1 : -1)

    return WaveState(modulation, amplitude, plasma, Ψ0, wave_direction, wave_along_b)
end

"""
    wave_frequency(modulation::AbstractWaveModulation, wave::WaveState, t)
    wave_frequency(wave::WaveState, t)

Calculates the frequency `ω` at some time `t` given some `WaveState`.

# Arguments

- `modulation`: Struct holding the modulation of the wave
- `wave`: Struct holding the defined state of the wave.
- `t`: Time in wave-duration to evaluate for [s].

# Returns

- The frequency `ω` at the given time `t` [rad/s].
"""
function wave_frequency(modulation::ChirpModulation, wave::WaveState, t)
    frequency_drift_rate = (modulation.ω1 - modulation.ω0) / modulation.t_duration
    return clamp(modulation.ω0 + frequency_drift_rate * t, modulation.ω0, modulation.ω1)
end

function wave_frequency(modulation::MonochromaticModulation, wave::WaveState, t)
    return modulation.ω
end

# TODO: Finish this
function wave_frequency(modulation::DemekhovModulation, wave::WaveState, t)
    ω0 = wave.plasma.ω_lb[1]

    # Find the wave-vector for the lowest density
    i = argmin(wave.plasma.Ω_e)
    k = dispersion_relation_whistler_branch(
        ω0,
        wave.θ,
        wave.plasma.ω_pe[i],
        wave.plasma.Ω_e[i]
    )

    v_perp = 1.0  #TODO: Figure out how to get this

    Bw = amplitude(wave.amplitude_model, 0.0, 0.0)
    frequency_drift_rate = (abs(qₑ) * k * v_perp * Bw) / (2π * mₑ)

    return clamp(ω0 + frequency_drift_rate * t, ω0, wave.plasma.ω_lb[end])
end

function wave_frequency(wave::WaveState, t)
    return wave_frequency(wave.modulation, wave, t)
end



"""
    wave_launch_time(modulation::AbstractWaveModulation, ω_grid)

Calculates the launching time for some frequency `ω`, given some modulation.

# Arguments

- `modulation`: Struct holding the modulation of the wave.
- `ω_grid`: Frequencies to evaluate for [rad/s].

# Returns

- The time it takes for each frequency to be launched [s].
"""
function wave_launch_time(modulation::MonochromaticModulation, ω_grid)
    # All frequencies launch at t = 0.0
    return zeros(length(ω_grid))
end

function wave_launch_time(modulation::ChirpModulation, ω_grid)
    chirp_rate = (modulation.ω1 - modulation.ω0) / modulation.t_duration
    iszero(chirp_rate) && return zeros(length(ω_grid))
    return [(ω - modulation.ω0) / chirp_rate for ω in ω_grid]
end

function wave_launch_time(modulation::DemekhovModulation, ω_grid)
    throw("Demekhov wave is not fully implemented yet...")
end

function wave_launch_time(wave::WaveState, ω_grid)
    return wave_launch_time(wave.modulation, ω_grid)
end

function wave_launch_time(wave::Nothing, ω_grid; t_duration=0.2)
    chirp_rate = (ω_grid[end] - ω_grid[1]) / t_duration
    iszero(chirp_rate) && return zeros(length(ω_grid))
    return [(ω - ω_grid[1]) / chirp_rate for ω in ω_grid]
end



"""
    wave_fields(
    wave::WaveState,
    t::Float64,
    s::Float64,
    Ψ::Float64,
    B0::SVector{3, Float64},
    λ::Float64
)

Calculates perturbations to the magnetic- and electric field, `δB` and `δE`.

The method is described by Saito et al. (2012). The wave is a right-hand polarized whistler
wave, and the cold plasma approximation is used to link the perturbations calculated in the
magnetic field to the electric field.

# Arguments

- `wave`: Struct holding the state of the wave.
- `t`: Current time [s].
- `s`: Current position along the field line, used for position-dependent amplitudes [m].
- `Ψ`: Current wave phase [rad].
- `B0`: Magnetic field at initial position, constructs the magnetic basis vectors [T].
- `λ`: Current magnetic latitude [rad].
"""
function wave_fields(
    wave::WaveState,
    t::Float64,
    s::Float64,
    Ψ::Float64,
    B0::SVector{3, Float64},
    λ::Float64
)

    # Get the frequency- and amplitude evolution of the wave based on the given models
    ω = wave_frequency(wave, t)
    Bw = wave_amplitude(wave.amplitude, t, s)

    _, e1, e2, _ = magnetic_basis(B0)

    # Find the latitude (and thus plasma parameters) in the λ_grid of the plasma closest to the one we're looking for
    i = searchsortedfirst(wave.plasma.λ_grid, λ)
    i = clamp(i, 1, length(wave.plasma.Ω_e))
    Ω_e = wave.plasma.Ω_e[i]
    ω_pe = wave.plasma.ω_pe[i]

    # Use the dispersion relation to find the SIGNED wave number of this position
    θ = 0.0     # Wave-normal angle, approximation is only valid for non-oblique waves.
    k_parallel = wave.wave_along_b * dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)

    # Find the field perturbations from equations 10 and 11 in Saito 2012
    # TODO: Check polarization and directions compared to the paper!
    δB = Bw * (cos(Ψ) .* e1 .- sin(Ψ) .* e2)

    Ew = (ω/k_parallel) * Bw   # From the cold plasma approximation
    δE = Ew * (-sin(Ψ) .* e1 .- cos(Ψ) .* e2)

    return (δB=δB, δE=δE, k_parallel=k_parallel, ω=ω)
end


function Base.show(io::IO, w::WaveState)
    println(io, "WaveState:")
    println(io, "  Wave propagating " * (w.wave_direction>0 ? "northward" : "southward"))
    println(io, "  Modulation = $(w.modulation)")
    println(io, "  Amplitude  = $(w.amplitude)")
    println(io, "  Ψ₀         = $(w.Ψ0)")
end
