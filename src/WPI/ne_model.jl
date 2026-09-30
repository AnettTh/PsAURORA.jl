using AURORA
using CairoMakie

# Bin centers
L_centers = [3.5, 3.9, 4.4, 4.9, 5.5, 6.2, 7.0, 7.8]

# Bin edges (midpoints between centers, plus outer edges)
L_edges = [
    3.3,   # 3.5 - 0.2
    3.7,   # midpoint 3.5-3.9
    4.15,  # midpoint 3.9-4.4
    4.65,  # midpoint 4.4-4.9
    5.2,   # midpoint 4.9-5.5
    5.85,  # midpoint 5.5-6.2
    6.6,   # midpoint 6.2-7.0
    7.4,   # midpoint 7.0-7.8
    8.2,   # 7.8 + 0.4
]

# Table values
#              3.5,   3.9,   4.4,   4.9,  5.5,  6.2   7.0, 7.8
n_e0_vals = [530. , 380. , 230. , 140. , 83. , 39. , 15. , 7.7]  # [cm⁻³]
α_vals    = [  0.2,   0.4,   0.8,   0.9,  0.8,  1.3,  2.1, 1.6]

# These are not used for anything, but a part of the model??
# L_α_vals  = [  8.1,   5.9,   4.8,   5.2,  6.4,  5.5,  4.8, 6.1]


##
# TODO: Make version that traces the field-line, if the magnetic model is not dipolar
function ne_denton(L, λ; SI::Bool=true)

    R = L .* RE .* cos.(λ).^2

    # Find the right bin
    i = max(1, searchsortedfirst(L_edges, L) - 1)
    i = clamp(i, 1, length(n_e0_vals))

    # Convert to SI-units
    if SI
        n_e0 = n_e0_vals[i] * 1e6   # [m⁻³]
    else
        n_e0 = n_e0_vals[i]
    end

    α = α_vals[i]

    return n_e0 .* (L.*RE ./ R).^α
end
