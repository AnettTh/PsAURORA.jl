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


function subtracted_bimaxwellian(
    v_parallel,
    v_perp,
    N::Float64,
    T_parallel::Float64,
    A::Float64,
    α_lc::Float64,
    Δ::Float64,
    r::Float64,
    M::Int
)

    a2_parallel = 2 * T_parallel * eV_in_J / mₑ
    a_parallel = sqrt(a2_parallel)
    a2_perp = a2_parallel * (A + 1)
    v2_parallel = v_parallel.^2
    v2_perp = v_perp.^2

    # Helper for f0
    function f0(v2_parallel, v2_perp)
        prefactor = N / (π^(3/2) * a2_perp * a_parallel)
        exp_parallel = @. exp(- (v2_parallel)/(a2_parallel))
        exp_perp = @. exp(- (v2_perp) / (a2_perp))

        return @. prefactor * exp_parallel * exp_perp
    end

    v_dj = @. (5 * range(-M, M) * a_parallel) / (M)
    aj_perp = v_dj .* atan(α_lc)
    aj_parallel = fill(r * (v_dj[2] - v_dj[1]), length(-M:M))
    Njs, _, _, _ = compute_Nj(M, a_parallel, sqrt(a2_perp), α_lc, Δ, N, r)

    f_hat = 0.0

    for j in eachindex(v_dj)
        ajpa = aj_parallel[j]
        ajpe = aj_perp[j]
        iszero(Njs[j]) && continue
        f_hat += Njs[j] / (π^(3/2) * ajpe^2 * ajpa) *
                 exp(-(v_parallel - v_dj[j])^2 / ajpa^2) *
                 exp(-v_perp^2 / ajpe^2)
    end

    return f0(v2_parallel, v2_perp) .- f_hat

end


function compute_Nj(M, a_par, a_perp, α_LC, Δ, N, r)

    # The range of f̂ to use
    js = -M:M

    # Define range, drift velocity and drift velocity spacing
    vdjs = [5j * a_par / M for j in js]
    Δvd = vdjs[2] - vdjs[1]

    # Computing a for each j in the range, with floor to avoid singularity at j=0
    a_par_js  = [r * abs(Δvd) for _ in js]
    a_perp_min = 1e-10
    a_perp_js = [max(abs(vdj) * tan(α_LC), a_perp_min) for vdj in vdjs]

    # f0 at each drift point
    f0_at_vdj = [N / (π^(3/2) * a_perp^2 * a_par) * exp(-vdj^2 / a_par^2)
                 for vdj in vdjs]

    # Using f̂(v_dj, 0) ≈ (1 - Δ)f₀(v_dj, 0) to solve for Nⱼ
    Njs = [(1 - Δ) * f0_at_vdj[i] * π^(3/2) * a_perp_js[i]^2 * a_par_js[i]
           for i in 1:2M+1]

    # j=0: vd0=0, a_perp_0=0 → no contribution
    Njs[M+1] = 0.0

    return Njs, vdjs, a_par_js, a_perp_js
end
