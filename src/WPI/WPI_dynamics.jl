using AURORA

# TODO: Use this for TOF? Test for arrays of Ω_e and ω_pe
"""
    dispersion_relation_whistler_branch(ω, θ; ω_pe=37.9e3, Ω_e=9.48e3)

Calculate the wavenumber the whistler-branch of the Appleton-Hartree equation for oblique
waves (Hsieh 2022).

The function uses a defined whistler angular frequency `ω` and a wave-normal angle `θ` to
solve the dispersion relation. The electron plasma frequency `ω_pe` and the electron
cyclotron frequency `Ω_e` are functions of position, while `ω` is a function of time,
meaning that the function supports ranges either as a funciton of position or time, but not
both.

# Arguments

- `ω`: Wave angular frequency for the whistler-mode chorus wave [Hz].
- `θ`: Wave normal angle of the propagating wave, non-zero value indicates oblique wave
  [rad].

# Keyword Arguments

- `ω_pe`: Electron plasma frequencie, default is 37.9 kHz (Hsieh 2022) [Hz].
- `Ω_e`: Electron cyclotron frequencie, default is 9.48 kHz (Hsieh 2022) [Hz].

# Returns

- The wave number that satisfies the dispersion relation.

# Throws

- TODO: Add these when what they should be are known
"""
function dispersion_relation_whistler_branch(ω, θ; ω_pe=37.9e3, Ω_e=9.48e3)

    abs(θ) > 1 && throw(ArgumentError("Are you sure you are using radians for the WNA?"))

    ω > abs(0.5 * Ω_e) && throw(ArgumentError(
        "You are not in the LBC-range, sure this is right?"
        ))

    # TODO: Add cold-plasma check, I think that involves Ω_e and ω_pe??

    # NOTE: Can remove the throws for array-check when the rest of the code is safe
    if ω isa AbstractArray && Ω_e isa AbstractArray
        length(ω) == length(Ω_e) || throw(DimensionMismatch(
            "ω and Ω_e must have the same length when both are arrays — " *
            "they represent different physical dimensions (frequency vs latitude)"
        ))
    end
    if ω isa AbstractArray && ω_pe isa AbstractArray
        length(ω) == length(ω_pe) || throw(DimensionMismatch(
            "ω and ω_pe must have the same length when both are arrays"
        ))
    end

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
