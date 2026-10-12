using AURORA

abstract type AbstractDensityModel end

struct ConstantDensity <: AbstractDensityModel
    n_e :: Float64
end

struct DentonDensity <: AbstractDensityModel
    R_max :: Float64
    ϕ_eq  :: Float64
end

# Constructors for Denton density model
function DentonDensity(magnetic_field, L, ϕ_eq)
    R_max = find_R_max(magnetic_field, Spherical(L, 0.0, ϕ_eq))
    return DentonDensity(R_max, ϕ_eq)
end

function DentonDensity(::typeof(dipole_field), L, ϕ_eq)
    R_max = Float64(L) * RE
    return DentonDensity(R_max, ϕ_eq)
end


# Evaluate the methods
function electron_density(model::ConstantDensity, L, λ)
    return model.n_e
end

function electron_density(model::DentonDensity, L, λ)
    return denton_density_model(Spherical(L, λ, model.ϕ_eq), model.R_max)
end

# Make both callable as functions
(model::ConstantDensity)(L, λ) = electron_density(model, L, λ)
(model::DentonDensity)(L, λ) = electron_density(model, L, λ)
