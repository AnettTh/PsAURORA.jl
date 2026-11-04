using Dates
using TsyganenkoModels
using StaticArrays

# Models for estimating the magnetic field.
# Implements AbstractMagneticField with DipoleMagneticField and TsyganenkoMagneticField.

"""
    AbstractMagneticField

Abstract type for magnetic field models.
"""
abstract type AbstractMagneticField end


"""
    DipoleMagneticField

Dipole magnetic field model, callable with either `Cartesian` or `Dipolar` coordinates.
"""
struct DipoleMagneticField <: AbstractMagneticField end

(model::DipoleMagneticField)(r::Cartesian) = dipole_field(r)
(model::DipoleMagneticField)(r::Dipolar) = dipole_field(r)


"""
    TsyganenkoMagneticField

Tsyganenko magnetic field model, callable with either `Cartesian` or `Dipolar`
coordinates, also containing free parameters describing the solar wind conditions.

# Fields

- `time`: Given as a string in ISO8601 format.
- `pdyn`: Solar wind dynamic pressure [nPa].
- `dst`: Disturbance Storm Time index [nT].
- `byimf`: IMF By component [nT].
- `bzimf`: IMF Bz component [nT].
- `ps`: Dipole tilt angle [rad].
- `ϕ`: Magnetic longitude [rad].
- `component`: String to chose which component of the field to return; either external
  `ext`, internal `int` or the total field `both`.
"""
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

function TsyganenkoMagneticField(
    ϕ;
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

(model::TsyganenkoMagneticField)(position::Cartesian) = tsyganenko_field(
    position;
    time      = model.time,
    pdyn      = model.pdyn,
    dst       = model.dst,
    byimf     = model.byimf,
    bzimf     = model.bzimf,
    ps        = model.ps,
    component = model.component,
)

# TODO: Try to remove this, as it is not valid when λ≠0
(model::TsyganenkoMagneticField)(position::Dipolar) = tsyganenko_field(
    position;
    time      = model.time,
    pdyn      = model.pdyn,
    dst       = model.dst,
    byimf     = model.byimf,
    bzimf     = model.bzimf,
    ps        = model.ps,
    component = model.component,
)

function Base.show(io::IO, ::DipoleMagneticField)
    println(io, "DipoleMagneticField()")
end

function Base.show(io::IO, f::TsyganenkoMagneticField)
    println(io, "TsyganenkoMagneticField() at ϕ_eq=$(round(rad2deg(f.ϕ), digits=1))°:\n" *
           "  $(round(longitude_to_MLT(f.ϕ, degrees=false), digits=1)) MLT\n" *
           "  dst=$(f.dst) \n" *
           "  pdyn=$(f.pdyn) \n" *
           "  byimf=$(f.byimf) \n" *
           "  bzimf=$(f.bzimf) \n"
           )
end
