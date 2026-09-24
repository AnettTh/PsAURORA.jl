using AURORA

abstract type AbstractTOF end

struct SimpleTOF <: AbstractTOF
    z_distance :: Float64       # Straight-line distance [m]
end

struct FieldlineTOF <: AbstractTOF
    r_source :: Float64         # Stopping radius for backwards propagation, i.e. source [m]
end

struct WPISimpleTOF <: AbstractTOF
    plasma :: PlasmaState
    n      :: Int
    θ      :: Float64
end

struct WPIFieldlineTOF <: AbstractTOF
    plasma :: PlasmaState
    n      :: Int
    θ      :: Float64
end


# No wave, straight line
function time_of_flight(particle::ParticleState, model::SimpleTOF)
    return model.z_distance / (abs(particle.μ) * particle.v)
end

# No wave, Boris mover along field line
function time_of_flight(particle::ParticleState, model::FieldlineTOF)
    v0 = get_v0_from_Eμ(
        particle.magnetic_field,
        particle.r0,
        particle.E_eV,
        particle.μ;
        towards_equator=true
    )
    result = boris_mover_TOF(particle.magnetic_field, particle.r0, v0, model.r_source)
    isnothing(result) && throw(ArgumentError("Particle did not precipitate."))
    return result.tof
end

# WPI, field-independent
function time_of_flight(particle::ParticleState, ω, model::WPISimpleTOF)
    return WPI_TOF(ω, particle, model.plasma;
                   field_dependent=false, n=model.n, θ=model.θ)
end

# WPI, field-dependent
function time_of_flight(particle::ParticleState, ω, model::WPIFieldlineTOF)
    return WPI_TOF(ω, particle, model.plasma;
                   field_dependent=true, n=model.n, θ=model.θ)
end


## Test that it works:
pa = ParticleState(3e4, -cos(deg2rad(3)), [6.5*RE, 0.0, 0.0], dipole_field; relativistic=true)
pl = PlasmaState(range(0.0, deg2rad(50), length=500), 1.8e7, dipole_field, 6.5)

tof = time_of_flight(pa, SimpleTOF(6.5*RE))
tof = time_of_flight(pa, FieldlineTOF(6.5*RE))
tof = time_of_flight(pa, 0.4*pl.Ω_e[1], WPISimpleTOF(pl, 1, 0.0))
tof = time_of_flight(pa, 0.4*pl.Ω_e[1], WPIFieldlineTOF(pl, 1, 0.0))
