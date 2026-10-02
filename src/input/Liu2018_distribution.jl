using AURORA


## IDEA: Add option for several subtractions
function subtracted_bimaxwellian_Liu(
    v_parallel,
    v_perp,
    N::Float64,
    T_parallel::Float64,
    A::Float64,
    Δ::Float64,
    β::Float64
)

    a2_parallel = 2 * T_parallel * eV_in_J / mₑ
    a2_perp = a2_parallel * (A + 1)
    v2_parallel = v_parallel.^2
    v2_perp = v_perp.^2

    prefactor = N / (π^(3/2) * a2_perp * sqrt(a2_parallel))

    exp_parallel = @. exp(- (v2_parallel)/(a2_parallel))

    exp_perp1 = @. exp(- (v2_perp) / (a2_perp))
    exp_perp2 = @. exp(- (v2_perp) / (a2_perp * β))

    frac = (1 - Δ) / (1 - β)

    return @. prefactor * exp_parallel * (Δ * exp_perp1 + frac * (exp_perp1 - exp_perp2))

end
