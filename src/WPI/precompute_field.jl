using AURORA
using JLD2
using Interpolations

"""
    get_magnetic_field(
    magnetic_field::AbstractMagneticField,
    L,
    ϕ_eq;
    cache_dir="src/WPI/data/field_cache",
    n_points=1000
)

Function to calculate the magnetic field along a fieldline, if it not already exists.

Takes some magnetic field and traces its fieldline from the initial position in the
equatorial plane, given by `L` and `ϕ_eq`, interpolated to a vector of length `n_points`
with equal spacing. The result is then saved with a name containint all the parameters used,
such that if there alredy exists such a file, the computation is skipped and the field
values from the field is returned instead.

# Arguments

- `magnetic_field`: Magnetic field function `f(x, y, z)`.
- `L`: L-shell.
- `ϕ_eq`: Equatorial longitude [rad].

# Keyword Arguments

- `cache_dir`: Path to where the data is stored, default is `src/WPI/data/field_cache`.
- `n_points`: Length of the resulting evaluated field, default is `1000`.

# Returns

- A callable  `f(s)` that returns `[Bx, By, Bz]` [T] at arc length `s` [m] along a field
  line where `s = 0` corresponds to one ionospheric footpoint, `s = s_max` the other.
"""
function get_magnetic_field(
    magnetic_field::AbstractMagneticField,
    L,
    ϕ_eq;
    cache_dir="src/WPI/data/field_cache",
    n_points=1000
)

    # Make unique name for the file
    # NOTE: Might also need to add more info, such as tsyg kwargs
    # NOTE: Trying to account for floating-point error by rounding of L, may need to look into better solution?
    model = magnetic_field isa TsyganenkoMagneticField ? "tsyg" : "dipole"
    filename = joinpath(
        cache_dir,
        "$(model)_L$(round(L, digits=4))_phi$(round(rad2deg(ϕ_eq), digits=1))_n$(n_points).jld2"
    )

    # Check if the file exists, else compute it
    if isfile(filename)
        @info "Loading precomputed $(model) field at L=$(round(L, digits=4))," *
        " ϕ=$(round(rad2deg(ϕ_eq), digits=1))° using n=$(n_points) points."
        @load filename s_grid Bx_vals By_vals Bz_vals
    else
        @info "Computing magnetic field along field line..."
        mkpath(cache_dir)

        # Compute using appropriate field-line tracing
        if magnetic_field isa TsyganenkoMagneticField

            # Trace the field-line
            _, trace = find_R_max(
                magnetic_field,
                L,
                0.0,
                ϕ_eq;
                ds=RE*0.001,
                store_trace=true
            )
            xs, ys, zs = trace.xs, trace.ys, trace.zs


            Bs = [magnetic_field(Cartesian(x, y, z)) for (x, y, z) in zip(xs, ys, zs)]

            Bx_fine = [B[1] for B in Bs]
            By_fine = [B[2] for B in Bs]
            Bz_fine = [B[3] for B in Bs]

            # Compute the arc length, not uniform spacing
            ds_vals = sqrt.(diff(xs).^2 .+ diff(ys).^2 .+ diff(zs).^2)
            s_fine  = [0.0; cumsum(ds_vals)]

            # Resample on defined grid with size n_points, uniform spacing
            interpolated_Bx = linear_interpolation(s_fine, Bx_fine)
            interpolated_By = linear_interpolation(s_fine, By_fine)
            interpolated_Bz = linear_interpolation(s_fine, Bz_fine)

            s_grid = collect(range(0.0, s_fine[end], length=n_points))
            Bx_vals = interpolated_Bx.(s_grid)
            By_vals = interpolated_By.(s_grid)
            Bz_vals = interpolated_Bz.(s_grid)


            @save filename s_grid Bx_vals By_vals Bz_vals
            @info "Saved to $(filename)"

        else
            # Define the latitude grid based on the dipole approximation
            λ_max = acos(sqrt((RE + z_ionosphere) / (L * RE)))
            λ_grid_fine = range(-λ_max, λ_max, length=10*n_points)

            # Compute positions along dipole field line
            xs = [L * RE * cos(λ)^3 for λ in λ_grid_fine]
            ys = zeros(length(λ_grid_fine))
            zs = [L * RE * cos(λ)^2 * sin(λ) for λ in λ_grid_fine]

            # Find the arc length
            ds_vals = sqrt.(diff(xs).^2 .+ diff(ys).^2 .+ diff(zs).^2)
            s_fine = [0.0; cumsum(ds_vals)]

            # Compute B analytically
            Bs = [magnetic_field(Cartesian(x, 0.0, z)) for (x, z) in zip(xs, zs)]
            Bx_fine = [B[1] for B in Bs]
            By_fine = [B[2] for B in Bs]
            Bz_fine = [B[3] for B in Bs]

            # Resample on defined grid with size defined by n_points
            interpolated_Bx = linear_interpolation(s_fine, Bx_fine)
            interpolated_By = linear_interpolation(s_fine, By_fine)
            interpolated_Bz = linear_interpolation(s_fine, Bz_fine)

            s_grid  = collect(range(0.0, s_fine[end], length=n_points))
            Bx_vals = interpolated_Bx.(s_grid)
            By_vals = interpolated_By.(s_grid)
            Bz_vals = interpolated_Bz.(s_grid)


            @save filename s_grid Bx_vals By_vals Bz_vals
            @info "Saved to $(filename)"
        end
    end

    # Build final interpolator
    final_interpolated_Bx = linear_interpolation(s_grid, Bx_vals; extrapolation_bc=Flat())
    final_interpolated_By = linear_interpolation(s_grid, By_vals; extrapolation_bc=Flat())
    final_interpolated_Bz = linear_interpolation(s_grid, Bz_vals; extrapolation_bc=Flat())

    return s -> [
        final_interpolated_Bx(s),
        final_interpolated_By(s),
        final_interpolated_Bz(s)
        ], s_grid
end
