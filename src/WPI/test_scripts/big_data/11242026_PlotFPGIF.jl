using AURORA
using JLD2
using ProgressMeter
using Dates
using CairoMakie
using Printf

#@load "src/WPI/data/11-16-2026/energy-time-position.jld2"
@load "src/WPI/data/11-25-2026/energy-time-position.jld2"

##
ref = LTPReference(params.station)

E_keV = collect(E_grid) ./ 1e3

east = fp_enu[1, :, :, :, :] ./ 1e3      # n_E × n_μ × n_ω × n_Ψ, in km
north = fp_enu[2, :, :, :, :] ./ 1e3

# Energy of every entry, same shape as east and north
E_all = repeat(E_keV, 1, size(east, 2), size(east, 3), size(east, 4))

# Flatten and keep only precipitated cases (NaN = no resonance or bounced)
precipitated = findall(isfinite, vec(east))
Δeast = vec(east)[precipitated]
Δnorth = vec(north)[precipitated]
E = vec(E_all)[precipitated]
TOF_all = vec(TOF)[precipitated]

@show length(precipitated), length(east)      # number of footpoints out of the total


## Pool in time
n_frames = 60
t_edges  = range(minimum(TOF_all), maximum(TOF_all); length = n_frames + 1)
t_frame  = collect(t_edges[2:end])               # time shown in each frame

## Observables that the plot reads from
pts    = Observable(Point2f[])
cols   = Observable(Float32[])
t_text = Observable("")

## Time frames: each frame shows everything that has arrived by the end of its bin
n_frames = 40
t_edges  = range(minimum(TOF_all), maximum(TOF_all); length = n_frames + 1)
t_frame  = collect(t_edges[2:end])               # time shown in each frame

## Observables that the plot reads from
pts    = Observable(Point2f[])
cols   = Observable(Float32[])
t_text = Observable("")

## Make figure
fig = Figure(size = (700, 350))
ax = Axis(fig[1, 1];
    xlabel = "East [km]",
    ylabel = "North [km]",
    title  = "Footpoints relative to $(params.station) at $(params.time)",
    aspect = DataAspect(),
)
sc = scatter!(ax, pts;
    color      = cols,
    colormap   = :turbo,
    colorrange = extrema(E_keV),
    markersize = 6,
)
Colorbar(fig[1, 2], sc; label = "E [keV]")
text!(ax, 0.02, 0.96; text = t_text, space = :relative, align = (:left, :top), fontsize = 14)

ylims!(ax, 4.0, 5.0)
xlims!(ax, 0.0, 2.5)

## Record
mkpath("src/WPI/test_scripts/big_data/figures")
record(fig, "src/WPI/test_scripts/big_data/figures/11262026FootpointsInTime.gif", 1:n_frames; framerate = 4) do k
    # Keeps everything that has arrived
    arrived = TOF_all .<= t_frame[k]
    pts[]    = Point2f.(Δeast[arrived], Δnorth[arrived])
    cols[]   = Float32.(E[arrived])
    t_text[] = @sprintf("t = %6.2f s", t_frame[k])
end
