abstract type AbstractDensityModel end

struct ConstantDensity <: AbstractDensityModel
    n_e :: Float
end

struct DentonDensity <: AbstractDensityModel

end


function electron_density(model::ConstantDensity, L, λ, magnetic_field)
    return model.n_e
end

function electron_density(model::DentonDensity, L, λ, magnetic_field)

end
