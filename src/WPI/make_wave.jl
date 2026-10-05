using AURORA

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
