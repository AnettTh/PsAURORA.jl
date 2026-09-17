using AURORA

struct ParticleState{F<:Function}
    E_eV::Float64               # Energy [eV]
    v::Float64                  # Speed [m/s]
    μ::Float64                  # Pitch-angle cosine
    α0::Float64                  # Pitch-angle [rad]
    r0::SVector{3, Float64}     # Initial position [m]
    B0::SVector{3, Float64}     # Magnetic field at initial position [T]
    b̂::SVector{3, Float64}      # Field unit-vector
    L::Float64                  # L-shell
    B_eq::SVector{3, Float64}   # Magnetic field at equator [T]
    α_eq::Float64               # Pitch-angle at equator
    α_lc::Float64               # Loss-cone angle based on initial position [rad]
    magnetic_field::F           # Magnetic field model
end

# TODO: Figure out how to better solve r0 and L, such that r0 can be anywhere and L is used for anything related to equator
function ParticleState(E_eV, μ, r0, magnetic_field)

    v = velocity_from_kinetic_energy(E_eV, mₑ)
    B0 = magnetic_field(r0...)
    b̂, _, _, _ = magnetic_basis(B0)
    r_mag = norm(r0)
    λ0 = asin(r0[3] / r_mag)
    L = r_mag / (RE * cos(λ0)^2)
    α_lc = losscone_angle(magnetic_field, r0)
    α0 = acos(μ)
    B_eq = BE / L^3
    α_eq = pitch_angle_at_z(α0, norm(B0), norm(B_eq))

    return ParticleState(
        E_eV,
        v,
        μ,
        α0,
        SVector{3}(r0),
        B0,
        SVector{3}(b̂),
        L,
        B_eq,
        α_eq,
        α_lc,
        magnetic_field
    )
end
