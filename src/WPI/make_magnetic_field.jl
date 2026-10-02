using Dates
using TsyganenkoModels

## Coordinate types
struct Cartesian
    x :: Float64
    y :: Float64
    z :: Float64
end

Base.iterate(r::Cartesian) = (r.x, 1)
Base.iterate(r::Cartesian, i) = i == 1 ? (r.y, 2) : i == 2 ? (r.z, 3) : nothing
Base.length(::Cartesian) = 3
Base.getindex(r::Cartesian, i::Int) = i == 1 ? r.x : i == 2 ? r.y : i == 3 ? r.z : throw(BoundsError(r, i))


import LinearAlgebra: norm
norm(r::Cartesian) = sqrt(r.x^2 + r.y^2 + r.z^2)

struct Spherical
    L :: Float64
    λ :: Float64
    ϕ :: Float64
end


# Abstract magnetic field type
abstract type AbstractMagneticField end

# Make the dipole model
struct DipoleMagneticField <: AbstractMagneticField end

# Call it using given coordinate system
(model::DipoleMagneticField)(r::Cartesian) = dipole_field(r.x, r.y, r.z)
(model::DipoleMagneticField)(r::Spherical) = dipole_field(r.L, r.λ)


# Make the tsyganenko model
struct TsyganenkoMagneticField <: AbstractMagneticField
    time      :: DateTime
    pdyn      :: Float64
    dst       :: Float64
    byimf     :: Float64
    bzimf     :: Float64
    ps        :: Float64
    ϕ         :: Float64
    component :: String
end

# Call it using given coordinate system
function TsyganenkoMagneticField(ϕ;
    time      = DateTime("2020-01-01T00:01:40"),
    pdyn      = 2.0,
    dst       = -87.0,
    byimf     = 2.0,
    bzimf     = -5.0,
    ps        = -0.533585131,
    component = "both"
)
    return TsyganenkoMagneticField(time, pdyn, dst, byimf, bzimf, ps, ϕ, component)
end

(model::TsyganenkoMagneticField)(r::Cartesian) = tsyganenko_field(r.x, r.y, r.z;
    time      = model.time,
    pdyn      = model.pdyn,
    dst       = model.dst,
    byimf     = model.byimf,
    bzimf     = model.bzimf,
    component = model.component
)
(model::TsyganenkoMagneticField)(r::Spherical) = tsyganenko_field_spherical(r.L, r.λ, r.ϕ;
    time      = model.time,
    pdyn      = model.pdyn,
    dst       = model.dst,
    byimf     = model.byimf,
    bzimf     = model.bzimf,
    component = model.component
)
