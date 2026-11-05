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

# -------------------------Table values------------------------
#              3.5,   3.9,   4.4,   4.9,  5.5,  6.2   7.0, 7.8
# -------------------------------------------------------------
n_e0_vals = [530. , 380. , 230. , 140. , 83. , 39. , 15. , 7.7]     # [cm⁻³]
α_vals    = [  0.2,   0.4,   0.8,   0.9,  0.8,  1.3,  2.1, 1.6]
L_α_vals  = [  8.1,   5.9,   4.8,   5.2,  6.4,  5.5,  4.8, 6.1]     # Scale length


"""
    denton_density_model(position::Cartesian, R_max, L_Denton_bin; SI::Bool=true)
    denton_density_model(position::Dipolar, R_max, L_Denton_bin; SI::Bool=true)

Compute density along a fieldline, given some initial position.

The model is implemented as described by Denton 2002. The initial position is given either
in spherical coordinates `L, λ, ϕ` or Cartesian coordinates `x, y, z` and a precomputed
maximum radial distance `R_max` is found by evaluating the magnetic field model. The density
is by default returned in SI units [m⁻³], but can also be given in #/cc [cm⁻³].

# Agurments

- `position`: Current position of type `Cartesian` or `Dipolar`.
- `R_max`: Maximal distance from the Earth along the current field line [m].
- `L_Denton_bin`: L-value deciding which bin to use, which depends on the chosen model
  (it will be the same as `R_max/RE` if the original model is used, but it will be decided
  by footpoints from fieldline tracing if the adapted model is to be used).

# Keyword Arguments

- `SI`: Return density in SI units or #/cc, default is `true`.

# Returns

- The density at the current position, using the given bin.

# Throws

- `ArgumentError`: If `L_Denton_bin` is out of the valid range for the model.
- `ArgumentError`: If λ is non-zero, as this requires other solutions when combined with the
  Tsyganenko field model.
"""
function denton_density_model(position::Cartesian, R_max, L_Denton_bin; SI::Bool=true)

    L_edges[1] ≤ L_Denton_bin ≤ L_edges[end] || throw(
        ArgumentError("L=$(round(L_Denton_bin, digits=2)) is outside of the valid range" *
        " [$(L_edges[1]), $(L_edges[end])]")
    )

    # Find current radius
    R  = @. sqrt(position.x^2 + position.y^2 + position.z^2)

    # Find correct bin
    i = clamp(searchsortedfirst(L_edges, L_Denton_bin) - 1, 1, length(n_e0_vals))

    # Store and return value based on the bin
    n_e0 = SI ? n_e0_vals[i] * 1e6 : n_e0_vals[i]
    α = α_vals[i]

    return n_e0 * (R_max / R)^α
end

function denton_density_model(position::Dipolar, R_max, L_Denton_bin; SI::Bool=true)
    # NOTE: Remove this when not needed anymore
    #abs(position.λ) > 1e-6 &&
    #    throw(ArgumentError("Issue using Dipolar coordinates for of-equatorial latitudes."))
    return denton_density_model(Cartesian(position), R_max, L_Denton_bin; SI=SI)
end
