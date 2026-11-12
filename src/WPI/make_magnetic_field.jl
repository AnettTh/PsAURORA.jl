using Dates
using TsyganenkoModels
using StaticArrays

# Models for estimating the magnetic field.
# Implements AbstractMagneticField with DipoleMagneticField and TsyganenkoMagneticField.
# TODO: Remove ϕ-call all together, in favor for t?
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
- `component`: String to chose which component of the field to return; either external
  `ext`, internal `int` or the total field `both`.
"""
struct TsyganenkoMagneticField <: AbstractMagneticField
    t         :: DateTime
    pdyn      :: Float64
    dst       :: Float64
    byimf     :: Float64
    bzimf     :: Float64
    ps        :: Float64
    component :: String
    r_limit   :: Float64
end


function TsyganenkoMagneticField(
    t::DateTime;
    pdyn = 2.0,
    dst = -50.0,
    byimf = 1.0,
    bzimf = 2.0,
    ps = nothing,
    component = "both",
    r_limit = 15*RE,
)
    ps_used = isnothing(ps) ? dipole_tilt(t) : Float64(ps)
    return TsyganenkoMagneticField(
        t,
        pdyn,
        dst,
        byimf,
        bzimf,
        ps_used,
        component,
        r_limit
    )
end

function Base.show(io::IO, ::DipoleMagneticField)
    println(io, "DipoleMagneticField()")
end

function Base.show(io::IO, f::TsyganenkoMagneticField)
    println(io, "TsyganenkoMagneticField() at $(f.t)\n" *
           "  dst   = $(f.dst) \n" *
           "  pdyn  = $(f.pdyn) \n" *
           "  byimf = $(f.byimf) \n" *
           "  bzimf = $(f.bzimf) \n" *
           "  ps    = $(round(f.ps, digits=9))"
           )
end
