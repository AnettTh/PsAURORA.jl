struct ParticleState
    E_eV::Float64           # Energy [eV]
    v::Float64              # Speed [m/s]
    μ::Float64              # Pitch-angle cosine
    α::Float64              # Pitch-angle [rad]
    r0::SVector{3, Float64} # Position [m]
    B::SVector{3, Float64}  # Magnetic field [T]
    b̂::SVector{3, Float64}  # Field unit-vector
    L::Float64              # L-shell
    α_lc::Float64           # Loss-cone angle [rad]
end


function ParticleState(E_eV, μ, r0, magnetic_field)
    v = velocity_from_kinetic_energy(E_eV, mₑ)
    B = magnetic_field(r0...)
    b̂, _, _, _ = magnetic_basis(B)
    r_mag = norm(r0)
    λ = asin(r0[3] / r_mag)
    L = r_mag / (RE * cos(λ)^2)
    α_lc = losscone_angle(magnetic_field, r0)
    α = acos(μ)
    return ParticleState(E_eV, v, μ, α, SVector{3}(r0), B, b̂, L, α_lc)
end
