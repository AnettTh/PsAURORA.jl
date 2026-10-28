using AURORA
using StaticArrays

# TODO: Finish and test this!
# Define different models for the modulation of the wave
abstract type AbstractWaveModulation end
abstract type AbstractAmplitude end


# Wave with constant frequency
struct Monochromatic <: AbstractWaveModulation
    ω                :: Float64    # Wave frequency [rad/s]
end

# Chirp wave with a linear increase/decrease in frequency
struct Chirp <: AbstractWaveModulation
    ω0         :: Float64   # Start frequency [rad/s]
    ω1         :: Float64   # End frequency [rad/s]
    t_duration :: Float64   # Duration of chirp [s]
end

# Chirp wave as described by Demekhov 2016
struct Demekhov <: AbstractWaveModulation end


# Define different models for the amplitude of the wave
struct ConstantAmplitude <: AbstractAmplitude
    Bw :: Float64
end

wave_amplitude(a::ConstantAmplitude, t, s) = a.Bw    # Is the same no matter the time t and the position s

# TODO: Add also Gaussian-like amplitude from Chen 2020


# Collect everything to define the wave state
struct WaveState{M<:AbstractWaveModulation, A<:AbstractAmplitude}
    modulation     :: M
    amplitude      :: A
    plasma         :: PlasmaState      # Used for dispersion
    Ψ0             :: Float64          # Initial wave phase [rad]
    wave_direction :: Int              # +1 for northward, -1 for southward
end

# NOTE: rand() here might be an issue?
function WaveState(modulation, amplitude, plasma; Ψ0=2π*rand(), wave_direction=-1)
    return WaveState(modulation, amplitude, plasma, Ψ0, wave_direction)
end

# TODO: Wave already store modulation, should this be changed?
function wave_frequency(modulation::Chirp, wave::WaveState, t)
    frequency_drift_rate = (modulation.ω1 - modulation.ω0) / modulation.t_duration
    return clamp(modulation.ω0 + frequency_drift_rate * t, modulation.ω0, modulation.ω1)
end

function wave_frequency(modulation::Monochromatic, wave::WaveState, t)
    return modulation.ω
end

# TODO: This is now frequency as a function of time at one location, need time as a function of frequency
function wave_frequency(modulation::Demekhov, wave::WaveState, t)
    ω0 = wave.plasma.ω_lb[1]

    # Find the wave-vector for the lowest density
    i = argmin(wave.plasma.Ω_e)
    k = dispersion_relation_whistler_branch(ω0, wave.θ, wave.plasma.ω_pe[i], wave.plasma.Ω_e[i])

    v_perp = 1.0  #TODO: figure out how to get this

    Bw = amplitude(wave.amplitude_model, 0.0, 0.0)
    frequency_drift_rate = (abs(qₑ) * k * v_perp * Bw) / (2π * mₑ)

    return clamp(ω0 + frequency_drift_rate * t, ω0, wave.plasma.ω_lb[end])
end

function wave_frequency(wave::WaveState, t)
    return wave_frequency(wave.modulation, wave, t)
end


# All frequencies launch at the same time
function wave_launch_time(modulation::Monochromatic, ω_grid)
    return zeros(length(ω_grid))
end

function wave_launch_time(modulation::Chirp, ω_grid)
    chirp_rate = (modulation.ω1 - modulation.ω0) / modulation.t_duration
    return [(ω - modulation.ω0) / chirp_rate for ω in ω_grid]
end

function wave_launch_time(modulation::Demekhov, ω_grid)
    throw("Demekhov wave is not fully implemented yet...")
end

function wave_launch_time(wave::WaveState, ω_grid)
    return wave_launch_time(wave.modulation, ω_grid)
end

# NOTE: Temporary solution to no wave
function wave_launch_time(wave::Nothing, ω_grid; t_duration=0.2)
    chirp_rate = (ω_grid[end] - ω_grid[1]) / t_duration
    return [(ω - ω_grid[1]) / chirp_rate for ω in ω_grid]
end


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

    b̂, e1, e2, B_mag = magnetic_basis(B0)

    # Find the latitude (and thus plasma parameters) in the λ_grid of the plasma closest to the one we're looking for
    i = searchsortedfirst(wave.plasma.λ, λ)
    i = clamp(i, 1, length(wave.plasma.Ω_e))
    Ω_e = wave.plasma.Ω_e[i]
    ω_pe = wave.plasma.ω_pe[i]

    # Use the dispersion relation to find the wave number of this position
    k = dispersion_relation_whistler_branch(ω, 0.0, ω_pe, Ω_e)

    # Find the field perturbations from equations 10 and 11 in Saito 2012
    # TODO: Check polarization and directions compared to the paper!
    δB = Bw * (cos(Ψ) .* e1 .- sin(Ψ) .* e2)

    Ew = (ω/k) * Bw   # From the cold plasma approximation
    δE = Ew * (-sin(Ψ) .* e1 .- cos(Ψ) .* e2)

    return δB, δE
end
