using AURORA
using CairoMakie
using JLD2
using Dates
using ProgressMeter

## Define configuration
name = "src/WPI/data/11-09-2026/footpoint.jld2"
use_boris = true
relativistic = true
WPI = true
field_dependent = true
ns = [1, 2]
duct_enhancement = 0.0
station = :sodankyla

## Define magnetic field
t = DateTime("2020-12-21T02:00:00")
ref = LTPReference(:sodankyla)
pdyn = 2.0
dst = -28.0
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

## Define the particle
μ0 = -cos(deg2rad(0.5))

# Find the L-shell and equatorial position of the fieldline belonging to the station at the
# requested time
trace_L_shell = find_L_shell(magnetic_field, ref, t)
L = trace_L_shell.L
position_eq = trace_L_shell.position_eq

ϕ_eq = magnetic_field.ϕ
r0 = Cartesian(Dipolar(L, 0.0, ϕ_eq))
E_grid = range(7e3, 40e3, length=10)

particles = Vector{ParticleState}(undef, length(E_grid))
@showprogress for (i, E) in enumerate(E_grid)
    particles[i] = ParticleState(
        E,
        μ0,
        r0,
        magnetic_field,
        ϕ_eq;
        relativistic=relativistic,
    )
end

## Define the plasma
λ_grid = range(-deg2rad(50), deg2rad(50), length=200)
lb_low, lb_high = 0.1, 0.5
density = DentonDensity(
    magnetic_field,
    L,
    ϕ_eq;
    duct_enhancement=duct_enhancement,
)

plasma = PlasmaState(λ_grid, density; lb_low=lb_low, lb_high=lb_high)

## Define the wave
ωs = [minimum(plasma.Ω_e) * 0.3]
ω_frac = 0.3
amplitude = ConstantAmplitude(100e-12)
modulation = Chirp(plasma.ω_lb[1], plasma.ω_lb[end], 0.2)

wave = WPI ? WaveState(modulation, amplitude, plasma; wave_direction=-1) : nothing

## Calculate TOF and footpoints
TOF = zeros(length(E_grid), length(ωs), length(ns))
FP = fill(NaN, 3, length(E_grid), length(ωs), length(ns))

_, line = find_R_max(magnetic_field, position_eq; store_trace = true)

progress = Progress(length(particles))
Threads.@threads for i in eachindex(particles)
    for (j, nn) in enumerate(ns)
        result = try
            WPI_TOF(
                ωs,
                particles[i],
                plasma;
                wave=wave,
                use_boris=use_boris,
                field_dependent=field_dependent,
                n=nn,
                r_start=position_eq,
                get_footpoint=true,
                line=line
            )
        catch e
            @warn "WPI_TOF failed" i n=nn exception=e
            continue
        end
        TOF[i, :, j]     = result.tof
        FP[:, i, :, j]  .= result.footpoint
    end

    next!(progress)
end

## Convert to East-North-Up coordinates with chosen station as reference
reference = LTPReference(station)
fp_enu = fill(NaN, size(FP))
for j in axes(FP, 4), iω in axes(FP, 3), iE in axes(FP, 2)
    fp_enu[:, iE, iω, j] .= FP_to_ENU(FP[:, iE, iω, j], reference, t)
end


E_keV  = collect(E_grid) ./ 1e3
ω_norm = collect(ωs) ./ minimum(plasma.Ω_e)           # ω/Ω_eq, for the legend
markers = [:circle, :rect, :diamond, :utriangle, :star5, :cross]

east  = fp_enu[1, :, :, :] ./ 1e3                         # n_E × n_ω, in km
north = fp_enu[2, :, :, :] ./ 1e3

##
fig = Figure(size = (800, 300))
ax  = Axis(fig[1, 1];
    xlabel = "East [km]",
    ylabel = "North [km]",
    title  = "Footpoints relative to $(station)",
    aspect = DataAspect(),
)

# One scatter per frequency, coloured by energy
markers = Dict(1 => :circle, 2 => :utriangle)
sc = nothing
for (j, nn) in enumerate(ns)
    sc = scatter!(ax, east[:, 1, j], north[:, 1, j];
        color      = E_keV,
        colormap   = :turbo,
        colorrange = extrema(E_keV),
        marker     = markers[nn],
        markersize = 15,
    )
end
Colorbar(fig[1, 2], sc; label = "E [keV]")

# Legend for the marker shapes
Legend(fig[1, 3],
    [MarkerElement(color = :gray, marker = markers[nn], markersize = 12) for nn in ns],
    ["n = $(nn)" for nn in ns],
)

#ylims!(-0.8, -0.4)
#save("src/WPI/test_scripts/TOF/figures/footpoints_enu.png", fig)
