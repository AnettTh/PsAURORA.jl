using AURORA
using CairoMakie

## Define method
duct_enhancement = 0.5
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T10:00:00")
pdyn = 2.0
dst = -28.0
byimf = 0.0
bzimf = 2.0

tsyg_field = TsyganenkoMagneticField(
    t;
    pdyn=pdyn,
    dst=dst,
    byimf=byimf,
    bzimf=bzimf
)
dipole_field = DipoleMagneticField()

## Define the stations fieldline
ref = LTPReference(station)
tsyg_fieldline = StationFieldLine(tsyg_field, ref, t)
dipole_fieldline = StationFieldLine(dipole_field, ref, t)

## Define the particles
α0 = 2.0
μ0 = -cos(deg2rad(α0))
E = 10e3

particle_tsyg = ParticleState(E, μ0, tsyg_fieldline)
particle_dipole = ParticleState(E, μ0, dipole_fieldline)

## Make a plasma
density_tsyg = DentonDensity(tsyg_fieldline)
density_dipole = DentonDensity(dipole_fieldline)

plasma_tsyg = PlasmaState(tsyg_fieldline, density_tsyg)
plasma_dipole = PlasmaState(dipole_fieldline, density_dipole)

## Frequencies
ω_frac = range(0.1, 0.5, length=100)
ω_tsyg = minimum(plasma_tsyg.Ω_e) .* ω_frac
ω_dipole = minimum(plasma_dipole.Ω_e) .* ω_frac


## Find minimal resonance energy
E_min_tsyg = minimal_resonance_energy(ω_tsyg, plasma_tsyg, particle_tsyg)
E_min_dipole = minimal_resonance_energy(ω_dipole, plasma_dipole, particle_dipole)

## Make the figure
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="ω [rad/s]",
    ylabel="Eₘᵢₙ [keV]",
    title="Threshold for cyclotron resonance at $station $t",
    #yscale=log10,
)

lines!(ax, ω_frac, E_min_tsyg / eV_in_J / 1e3, label="Tsyganenko")
lines!(ax, ω_frac, E_min_dipole / eV_in_J / 1e3, label="Dipole")

axislegend(ax; position=:rt)
xlims!(ax, 0.1, 0.5)
ylims!(ax, 7, 40)       # Realistic values

save("src/WPI/test_scripts/minimal_resonance/figures/11192026_MinimalResonanceEnergy_$t.png", fig)
