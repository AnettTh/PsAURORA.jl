using AURORA
using StaticArrays

# TODO: Finish and test this!
abstract type AbstractWave end

struct MonochromaticWave <: AbstractWave
    ω                :: Float64    # Wave frequency [rad/s]
    θ                :: Float64    # Wave normal angle [rad]
    A_Bw             :: Float64    # Wave magnetic amplitude [T]
    electric_w_field :: Any        # PLACEHOLDER
    magnetic_w_field :: Any        # PLACEHOLDER
end

function MonochromaticWave(ω, θ, Bw=0.0)
    return MonochromaticWave(
        ω,
        θ,
        Bw,
        (x, y, z) -> zeros(3),
        (x, y, z) -> zeros(3)
    )
end

struct ChirpWave <: AbstractWave
    ω0         :: Float64   # Start frequency [rad/s]
    ω1         :: Float64   # End frequency [rad/s]
    t_duration :: Float64   # Duration of chirp [s]
    θ          :: Float64   # Wave normal angle [rad]
    A_Bw       :: Float64   # Wave magnetic amplitude [T]
    E_field    :: Any       # PLACEHOLDER
    B_field    :: Any       # PLACEHOLDER
end

function ChirpWave(ω0, ω1, t_duration, θ, A_Bw)
    return ChirpWave(ω0, ω1, t_duration, θ, A_Bw,
        (x, y, z) -> zeros(3),      # PLACEHOLDER
        (x, y, z) -> zeros(3))      # PLACEHOLDER
end


struct DemekhovWave <: AbstractWave
    plasma     :: PlasmaState   # Holding the plasma parameters
    t_duration :: Float64       # Duration of chirp [s]
    θ          :: Float64       # Wave normal angle [rad]
    A_Bw       :: Float64       # Wave magnetic amplitude [T]
    E_field    :: Any           # PLACEHOLDER
    B_field    :: Any           # PLACEHOLDER
end

function DemekhovWave(plasma::PlasmaState, t_duration, θ, A_Bw)
    return DemekhovWave(plasma, t_duration, θ, A_Bw,
        (x, y, z) -> zeros(3),      # PLACEHOLDER
        (x, y, z) -> zeros(3))      # PLACEHOLDER
end

function wave_frequency(wave::MonochromaticWave, t)
    frequency_drift_rate = 0
    return wave.ω + frequency_drift_rate * t
end

function wave_frequency(wave::ChirpWave, t)
    frequency_drift_rate = (wave.ω1 - wave.ω0) / wave.t_duration
    return clamp(wave.ω0 + frequency_drift_rate * t, wave.ω0, wave.ω1)
end

function wave_frequency(wave::DemekhovWave, t)
    # TODO: This is now frequency as a function of time at one location, need time as a function of frequency
    ω0 = wave.plasma.ω_lb[1]
    k = dispersion_relation_whistler_branch(ω0, wave.θ, wave.plasma.ω_pe[1], wave.plasma.Ω_e[1])
    v_g = group_velocity_whistler_wave(ω0, wave.plasma.Ω_e[1], wave.plasma.ω_pe[1])

    frequency_drift_rate = (qₑ * k * v_g * wave.A_Bw) / (2π * mₑ)

    return clamp(ω0 + frequency_drift_rate * t, ω0, wave.plasma.ω_lb[end])
end


abstract type AbstractAmplitude end

struct ConstantAmplitude <: AbstractAmplitude
    Bw :: Float64
end

# TODO: Add also Gaussian-like amplitude from Chen 2020

amplitude(a::ConstantAmplitude, t, s) = a.Bw    # Is the same no matter the time t and the position s

# TODO: Dispatch also Gaussian-like amplitude when implemented



function wave_fields(
    wave::AbstractWave,
    amplitude_model::AbstractAmplitude,
    t::Float64,
    s::Float64,
    Ψ::Float64,
    B0::SVector{3, Float64},
    plasma::PlasmaState,
    λ::Float64
)

    # Get the frequency- and amplitude evolution of the wave based on the given models
    ω = wave_frequency(wave, t)
    Bw = amplitude(amplitude_model, t, s)

    b̂, e1, e2, B_mag = magnetic_basis(B0)

    # Find the latitude (and thus plasma parameters) in the λ_grid of the plasma closest to the one we're looking for
    i = searchsortedfirst(plasma.λ, λ)
    i = clamp(i, 1, length(plasma.Ω_e))
    Ω_e = plasma.Ω_e[i]
    ω_pe = plasma.ω_pe[i]

    # Use the dispersion relation to find the wave number of this position
    k = dispersion_relation_whistler_branch(ω, 0.0, ω_pe, Ω_e)

    # Find the field perturbations from equations 10 and 11 in Saito 2012
    # TODO: Check polarization and directions compared to the paper!
    δB = Bw * (cos(Ψ) .* e1 .- sin(Ψ) .* e2)

    Ew = (ω/k) * Bw   # From the cold plasma approximation
    δE = Ew * (-sin(Ψ) .* e1 .- cos(Ψ) .* e2)

    return δB, δE
end
