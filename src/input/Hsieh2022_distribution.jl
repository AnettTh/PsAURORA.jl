using AURORA

# TODO: This also has opposite behavior, figure out why
function lorentz_factor(V_r, α_h)
    return c₀ / sqrt(c₀^2 - (V_r^2 * (sec(α_h))^2))
end

"""
    subtracted_bimaxwellian_Hsieh(
    V_r::Real,
    α_h::Real,
    v_parallel,
    v_perp,
    n_h:Real,
    V_parallel::Real,
    V_perp::Real,
    ρ::Real,
    β::Real
)

TBW
"""
function subtracted_bimaxwellian_Hsieh(
    K_J,
    α,
    n_h::Real,
    V_t_parallel::Real,
    V_t_perp::Real,
    ρ::Real,
    β::Real
)
    E₀ = mₑ * c₀^2
    u2 = c₀^2 * ((K_J^2 / E₀^2) + (2*K_J / E₀))
    u = sqrt(u2)

    γ_t_parallel  = 1 / sqrt(1 - V_t_parallel^2  / c₀^2)
    γ_t_perp = 1 / sqrt(1 - V_t_perp^2 / c₀^2)
    U_t_parallel  = γ_t_parallel  * V_t_parallel
    U_t_perp = γ_t_perp * V_t_perp

    prefactor = (sqrt(2/π) * n_h * c₀^3) / (U_t_parallel * U_t_perp^2 * (1 - ρ * β))

    exp_parallel = @. exp(-(u2 * (cos(α))^2) / (2 * U_t_parallel^2))

    exp1_perp = @. exp(-(u2 * (sin(α))^2)/(2 * U_t_perp^2))
    exp2_perp = @. exp(-(u2 * (sin(α))^2)/(2 * β * U_t_perp^2))

    bracket = @. exp1_perp - (ρ * exp2_perp)

    jacobian = @. (u / c₀) * ((K_J / E₀^2) + (1 / E₀) ) * sin(α)

    return @. prefactor * exp_parallel * bracket * jacobian
end
