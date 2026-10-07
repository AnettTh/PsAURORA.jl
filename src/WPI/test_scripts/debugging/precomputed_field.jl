using AURORA
using LinearAlgebra
using JLD2

## Make magnetic field
ϕ_eq = deg2rad(240.0)
magnetic_field = TsyganenkoMagneticField(ϕ_eq)

## Make particle
L = 6.5
ϕ = deg2rad(240.0)
μ = -cos(deg2rad(2.0))
E = 30e3
r0 = Cartesian(Spherical(L, 0.0, ϕ))

p = ParticleState(
    E,
    μ,
    r0,
    magnetic_field,
    ϕ_eq;
    precomputed_magnetic_field=true,
    relativistic=true,
    φ0=deg2rad(0.0)
    )


## Test the particle
@show p.E_eV / 1e3           # should be 30 keV
@show rad2deg(p.α0)          # should be 2°
@show rad2deg(p.α_eq)        # should be ≈ α0 since we start at equator
@show rad2deg(p.α_lc)        # should be ~3-5° for L=6.5
@show p.α_eq < p.α_lc        # should be FALSE — particle outside loss cone
@show p.L                    # should be 6.5
@show p.B_eq * 1e9           # equatorial B in nT, should be ~100 nT for L=6.5
@show norm(p.B0) * 1e9       # field at r0 in nT, should ≈ B_eq
@show p.γ                    # Lorentz factor, ~1.06 for 30 keV
@show norm(p.v0) / p.v       # should be 1.0


## Test the interpolated magnetic field
B_direct = norm(magnetic_field(Spherical(L, 0.0, ϕ)))
B_interp = norm(p.B_interpolated(p.s0))
@show B_direct * 1e9   # nT
@show B_interp * 1e9   # should match B_direct
@show abs(B_direct - B_interp) / B_direct   # relative error, should be < 1%


## Test initial position at arc
s_grid  = range(0, 1, length=100) .* (2 * p.s0)   # rough range
B_along = [norm(p.B_interpolated(s)) for s in s_grid]
s_min   = s_grid[argmin(B_along)]
@show p.s0 / RE        # arc length at equator in RE
@show s_min / RE       # should match p.s0


## Test v0
b̂ = p.b̂
v_par  = dot(p.v0, b̂)          # parallel component
v_perp = norm(p.v0 - v_par*b̂)  # perpendicular component
@show v_par / p.v               # should be ≈ cos(α0) ≈ cos(3°) ≈ 0.999
@show v_perp / p.v              # should be ≈ sin(α0) ≈ sin(3°) ≈ 0.052


## Check minimum B in the interpolator
B_itp, s_grid = get_magnetic_field(magnetic_field, round(L, digits=4), ϕ)

B_along = [norm(B_itp(s)) for s in s_grid]
s_min   = s_grid[argmin(B_along)]

@show p.s0 / RE         # where ParticleState thinks equator is
@show s_min / RE        # where interpolator minimum is
@show minimum(B_along) * 1e9   # minimum B in interpolator
@show maximum(B_along) * 1e9   # maximum B in interpolator
@show B_along[1] * 1e9         # B at s=0
@show B_along[end] * 1e9       # B at s=s_max


##
# Where does the trace start?
B_itp, s_grid = get_magnetic_field(magnetic_field, round(L, digits=4), ϕ)

# Load trace positions
@load "src/WPI/data/field_cache/tsyg_L6.5_phi240.0_n1000.jld2" s_grid Bx_vals By_vals Bz_vals

# The equatorial point in the trace
idx_eq = argmin(sqrt.(Bx_vals.^2 .+ By_vals.^2 .+ Bz_vals.^2))
@show idx_eq
@show s_grid[idx_eq] / RE
