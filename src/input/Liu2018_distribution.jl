using AURORA

# TODO: This has opposite behavior, figure out why
##
function subtracted_bimaxwellian_Liu(
    v_parallel,
    v_perp,
    N::Float64,
    a_parallel::Float64,
    a_perp::Float64,
    Δ::Float64,
    β::Float64;
    relativistic::Bool=true
)

    if relativistic
        prefactor = 1

        parallel = 1

        perpendicular = 1
    else
        prefactor = N / ((π)^(3/2) * a_perp^2 * a_parallel)

        parallel = @. exp(-(v_parallel^2)./(a_parallel^2))

        av = @. (v_perp.^2)/(a_perp^2)
        perpendicular = @. Δ * exp(-av) + ((1 - Δ)/(1 - β)) * (exp(-av) - exp(- av / β))
    end

    return prefactor * parallel .* perpendicular
end
