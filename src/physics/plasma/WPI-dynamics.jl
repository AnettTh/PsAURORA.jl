using AURORA

# TODO: Add time-of-flight considering also WPI
# TODO: Add pitch-angle scattering as well
function electron_cyclotron_frequency(r; magnetic_field::Function=dipole_field)

    B = magnetic_field(r...)

    Ω_e = norm(eV_in_J * B / mₑ)

    return Ω_e
end


# TODO: Find some nice default plasma-parameters
function group_velocity_whistler_wave(ω; ω_pe=1, Ω_e=1)

    a = (2 * c) / (ω_pe / Ω_e)
    b = (1 - (ω / Ω_e))^(3/2)
    c = (ω / Ω_e)^(1/2)

    return a*b*c
end


# TODO: Make this into a look-up table in some clever way??
function B_of_z(z, magnetic_field)

end


function kz_of_z(z)
end


function resonance_energy(z, ω, B_of_z, kz_of_z; θ=0)

    kz = kz_of_z(ω, θ)

end




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
