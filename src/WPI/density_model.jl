using AURORA

## TODO: Look into returning  magnetic_field and ϕ also from this, to avoid pushing it twice!
abstract type AbstractDensityModel end

struct ConstantDensity <: AbstractDensityModel
    n_e :: Float64
end

struct DentonDensity <: AbstractDensityModel
    R_max :: Float64
    ϕ     :: Float64
end

# Constructors for Denton density model
function DentonDensity(magnetic_field, L, ϕ)
    R_max = find_R_max(magnetic_field, Float64(L), 0.0, ϕ)      # NOTE: should this be ϕ_eq?
    return DentonDensity(R_max, ϕ)
end

function DentonDensity(::typeof(dipole_field), L, ϕ)
    R_max = Float64(L) * RE
    return DentonDensity(R_max, ϕ)
end


# Evaluate the methods
function electron_density(model::ConstantDensity, L, λ)
    return model.n_e
end

function electron_density(model::DentonDensity, L, λ)
    return denton_density_model(L, λ, model.ϕ, model.R_max)
end

# Make both callable as functions
(model::ConstantDensity)(L, λ) = electron_density(model, L, λ)
(model::DentonDensity)(L, λ) = electron_density(model, L, λ)
