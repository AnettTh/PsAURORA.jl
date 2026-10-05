using AURORA

"""
    dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)

Calculate the wavenumber the whistler-branch of the Appleton-Hartree equation for oblique
waves (Hsieh 2022, equation ).

The function uses a defined whistler angular frequency `ω` and a wave-normal angle (WNA) `θ`
to solve the dispersion relation. The electron plasma frequency `ω_pe` and the electron
cyclotron frequency `Ω_e` are functions of position, while `ω` is a function of time,
meaning that the function supports ranges either as a funciton of position or time, but not
both.

# Arguments

- `ω`: Wave angular frequency for the whistler-mode chorus wave [rad/s].
- `θ`: Wave normal angle of the propagating wave, non-zero value indicates oblique wave
  [rad].
- `ω_pe`: Electron plasma frequencie [rad/s].
- `Ω_e`: Electron cyclotron frequencie [rad/s].

# Returns

- The wave number that satisfies the dispersion relation.

# Throws

- `ArgumentError`: If the WNA is outside of ±1 (not in radians).
- `ArgumentError`: If the frequency is outside the LBC range.
"""
function dispersion_relation_whistler_branch(ω, θ, ω_pe, Ω_e)

    abs(θ) > 1 && throw(ArgumentError("Are you sure you are using radians for the WNA?"))

    any(ω .> abs.(0.5 .* Ω_e)) && throw(ArgumentError(
        "You are not in the LBC-range, sure this is right?"
        ))

    X = @. ω_pe^2 / ω^2
    Y = @. Ω_e / ω

    a = @. 2*(1-X)
    sin2 = sin(θ)^2

    num = @. X*a
    denum = @. a - Y^2 * sin2 + Y*sqrt(Y^2 * sin2^2 + a^2 * cos(θ)^2)

    c2k2ω2 = @. 1 - (num/denum)

    k = @. sqrt(c2k2ω2) * ω / c₀

    return k
end

# NOTE: Might need a WNA here?
"""
    group_velocity_whistler_wave(ω, Ω_e, ω_pe)

Calculates the group velocity `v_g` for the whistler mode chorus wave as a function of
chorus angular frequency `ω`.

This is valid assuming chorus frequencies `ω` ≫ ion gyro-frequencies (Chen 2020).

# Arguments

- `ω`: Chorus angular frequency, often a linearly rising tone [rad/s].
- `Ω_e`: Electron gyrofrequency at latitude λ [rad/s].
- `ω_pe`: Electron plasma frequency at latitude λ [rad/s].

# Returns

- The group velocity for the set of parameters chosen [m/s].

# Throws

- `ArgumentError`: If the frequencies does not match that of the whistler-branch.
"""
function group_velocity_whistler_wave(ω, Ω_e, ω_pe)

    any(@. ω > Ω_e) && throw(ArgumentError("Not on whistler branch, check your frequencies!"))

    a = @. (2 * c₀) / (ω_pe / Ω_e)
    b = @. (1 - (ω / Ω_e))^(3/2)
    c = @. (ω / Ω_e)^(1/2)

    return @. a*b*c
end
