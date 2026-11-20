using AURORA
using JLD2
using ProgressMeter
using Dates


## Configurations; rest is default
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t     = DateTime("2020-12-21T02:00:00")
pdyn  = 2.0
dst   = -28.0
byimf = 0.0
bzimf = 2.0

magnetic_field = TsyganenkoMagneticField(t; pdyn=pdyn, dst=dst, byimf=byimf, bzimf=bzimf)

## Grids with TOF(E, μ0, ω, Ψ0) = (250, 5, 3, 4)
E_grid = range(7e3, 40e3, length=20)
α_grid = [0.5]
μ0_grid = -cos.(deg2rad.(α_grid))
ω_frac = [0.1, 0.3]
Ψ0_grid = range(0, 3π/2, step=π/2)

## Output
TOF = fill(NaN, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))
FP = fill(NaN, 3, length(E_grid), length(μ0_grid), length(ω_frac), length(Ψ0_grid))


## Define the stations fieldline
ref = LTPReference(:sodankyla)
fieldline = StationFieldLine(magnetic_field, ref, t)

## Make the plasma
density = DentonDensity(fieldline)
plasma = PlasmaState(fieldline, density)

## Define the frequencies sampeled
ωs = ω_frac .* minimum(plasma.Ω_e)

## Make the wave modulation for the full lower-band whistler mode chorus wave
modulation = ChirpModulation(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)
amplitude = ConstantAmplitude(100e-12)


## Make all combinations of energies and pitch-angles
index_pairs = CartesianIndices((length(E_grid), length(μ0_grid)))

## Build the grid of particles for each combination of E and μ
particles = Matrix{Union{Nothing, ParticleState}}(nothing, size(index_pairs)...)
progress  = Progress(length(index_pairs); desc="Building particles: ")

Threads.@threads :greedy for k in eachindex(index_pairs)
    iE, iμ = Tuple(index_pairs[k])
    particles[iE, iμ] = ParticleState(
        E_grid[iE],
        μ0_grid[iμ],
        fieldline;
    )
end

## Running the simulation for each wave- and particle combination
for (iΨ, Ψ0) in enumerate(Ψ0_grid)
    # Make wave with current Ψ0
    wave = WaveState(modulation, amplitude, plasma, fieldline; Ψ0=Ψ0)

    # Track the progress
    progress_run = Progress(length(index_pairs); desc="Ψ₀ = $(round(rad2deg(Ψ0), digits=0))°: ")

    # Run the TOF calculation for each of the energy- and pitch-angle pairs
    Threads.@threads :greedy for k in eachindex(index_pairs)
        iE, iμ = Tuple(index_pairs[k])
        p = particles[iE, iμ]
        result = WPI_TOF(
            ωs,
            p,
            fieldline,
            plasma;
            wave=wave,
            get_footpoint=true,
        )

        TOF[iE, iμ, :, iΨ] = result.tof
        FP[:, iE, iμ, :, iΨ] = result.footpoint

        next!(progress_run)
    end
end

##
ref = LTPReference(station)
fp_enu = fill(NaN, size(FP))

for iΨ in Base.axes(FP, 5), iω in Base.axes(FP, 4), iμ in Base.axes(FP, 3), iE in Base.axes(FP, 2)
    fp_enu[:, iE, iμ, iω, iΨ] .= FP_to_ENU(FP[:, iE, iμ, iω, iΨ], ref, t)
end

E_keV = collect(E_grid) ./ 1e3

east  = fp_enu[1, :, :, :, :] ./ 1e3      # n_E × n_μ × n_ω × n_Ψ, in km
north = fp_enu[2, :, :, :, :] ./ 1e3

# Energy of every entry, same shape as east and north
E_all = repeat(E_keV, 1, size(east, 2), size(east, 3), size(east, 4))

# Flatten and keep only precipitated cases (NaN = no resonance or bounced)
percipitated = findall(isfinite, vec(east))
Δeast = vec(east)[percipitated]
Δnorth = vec(north)[percipitated]
E = vec(E_all)[percipitated]

@show length(percipitated), length(east)      # number of footpoints out of the total

## Plot footpoints only
fig_FP = Figure(size = (700, 350))
ax_FP = Axis(fig_FP[1, 1];
    xlabel = "East [km]",
    ylabel = "North [km]",
    title  = "Footpoints relative to $(station) at $t \n" *
       "$(length(percipitated))/$(length(east)) precipitated",
    aspect = DataAspect(),
)

sc_FP = scatter!(
    ax_FP,
    Δeast,
    Δnorth;
    color      = E,
    colormap   = :turbo,
    colorrange = extrema(E_keV),
    markersize = 6,
)
Colorbar(fig_FP[1, 2], sc_FP; label = "E [keV]")
ylims!(4, 5)
xlims!(0, 2.5)
display(fig_FP)
