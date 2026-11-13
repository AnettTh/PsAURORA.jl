using AURORA
using JLD2
using ProgressMeter
using Dates
using WGLMakie
using Bonito

@load "src/WPI/data/11-13-2026/energy-time-position.jld2"

ref = LTPReference(params.station)
fp_enu = fill(NaN, size(FP))

for iΨ in Base.axes(FP, 5), iω in Base.axes(FP, 4), iμ in Base.axes(FP, 3), iE in Base.axes(FP, 2)
    fp_enu[:, iE, iμ, iω, iΨ] .= FP_to_ENU(FP[:, iE, iμ, iω, iΨ], ref, DateTime(params.time))
end

E_keV = collect(E_grid) ./ 1e3

east  = fp_enu[1, :, :, :, :] ./ 1e3
## Same shape as east: n_E × n_μ × n_ω × n_Ψ
tof   = TOF                                                  # [s]

valid = isfinite.(east) .& (tof .> 0)                        # precipitated cases only

# Energy of every entry, same shape as east
E_all = repeat(E_keV, 1, size(east, 2), size(east, 3), size(east, 4))

α_deg  = rad2deg.(collect(α_grid))
colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

##
fig = Figure(size = (900, 750))
ax = Axis3(fig[1, 1];
    xlabel = "Time of flight [s]",
    ylabel = "East [km]",
    zlabel = "E [keV]",
    title  = "Energy, time of flight and east position of the footpoints",
)

# One scatter per pitch angle: one colour and one legend entry each
for iμ in 1:size(east, 2)
    m = valid[:, iμ, :, :]                                   # mask for this pitch angle
    scatter!(ax,
        tof[:, iμ, :, :][m],                                 # x: time
        east[:, iμ, :, :][m],                                # y: east
        E_all[:, iμ, :, :][m];                               # z: energy
        color      = colors[mod1(iμ, length(colors))],
        markersize = 6,
        label      = "α = $(round(α_deg[iμ], digits=2))°",
    )
end

axislegend(ax; position = :rt)

# Serve the finished figure in the browser
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
