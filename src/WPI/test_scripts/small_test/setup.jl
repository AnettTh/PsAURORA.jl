using AURORA

use_boris        = true
field_dependent  = true
relativistic     = true
WPI              = true
adjusted_model   = false
n                = 1
duct_enhancement = 0.0
station          = :sodankyla
wave_direction   = -1

## Define magnetic field
t     = DateTime("2020-12-21T02:00:00")
ref   = LTPReference(:sodankyla)
pdyn  = 2.0
dst   = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(
    ref;
    time=t,
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)

##
