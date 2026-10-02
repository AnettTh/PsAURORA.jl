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
    ω0      :: Float64   # Start frequency [rad/s]
    ω1      :: Float64   # End frequency [rad/s]
    t       :: Float64   # Duration of chirp [s]
    θ       :: Float64   # Wave normal angle [rad]
    Bw      :: Float64   # Wave magnetic amplitude [T]
    E_field :: Any       # PLACEHOLDER
    B_field :: Any       # PLACEHOLDER
end

function ChirpWave(ω0, ω1, t, θ, Bw)
    return ChirpWave(ω0, ω1, t, θ, Bw, (x, y, z) -> zeros(3), (x, y, z) -> zeros(3))
end

function wave_frequency(wave::ChirpWave, t)
    chirp_rate = (wave.ω1 - wave.ω1) / wave.t
    return wave.ω0 + chirp_rate * t
end
