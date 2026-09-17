using AURORA
include("plasma_parameters.jl")



function dispersion_relation_whistler_branch(ω, θ; ω_pe=37.9e3, Ω_e=9.48e3)

    X = ω_pe^2 / ω^2
    Y = Ω_e / ω

    a = 2*(1-X)
    sin2 = sin(θ)^2

    num = X*a
    denum = a - Y^2 * sin2 + Y*sqrt(Y^2 * sin2^2 + a^2 * cos(θ)^2)

    c2k2ω2 = 1 - (num/denum)

    k = sqrt(c2k2ω2) * ω / c₀

    return k
end
