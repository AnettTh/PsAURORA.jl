using Dates
using TsyganenkoModels
using StaticArrays

# Abstract magnetic field type
abstract type AbstractMagneticField end

# Make the dipole model
struct DipoleMagneticField <: AbstractMagneticField end

# Call it using given coordinate system
(model::DipoleMagneticField)(r::Cartesian) = dipole_field(r.x, r.y, r.z)
(model::DipoleMagneticField)(r::Spherical) = dipole_field(Cartesian(r))


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
    dst       = -50.0,
    byimf     = 1.0,
    bzimf     = 2.0,
    ps        = -0.533585131,
    component = "both"
)
    return TsyganenkoMagneticField(time, pdyn, dst, byimf, bzimf, ps, ϕ, component)
end

(model::TsyganenkoMagneticField)(position::Cartesian) = tsyganenko_field(position;
    time      = model.time,
    pdyn      = model.pdyn,
    dst       = model.dst,
    byimf     = model.byimf,
    bzimf     = model.bzimf,
    ps        = model.ps,
    component = model.component,
)
(model::TsyganenkoMagneticField)(position::Spherical) = tsyganenko_field(position;
    time      = model.time,
    pdyn      = model.pdyn,
    dst       = model.dst,
    byimf     = model.byimf,
    bzimf     = model.bzimf,
    ps        = model.ps,
    component = model.component,
)
