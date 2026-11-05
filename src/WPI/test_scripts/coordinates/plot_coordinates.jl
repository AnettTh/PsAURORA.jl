using AURORA
using WGLMakie
using Bonito
using LinearAlgebra
using Dates

##
t = DateTime("2020-01-01T10:00:00")
x̂, ŷ, ẑ = sm_axes(t)
s = sun_unit_vector_GEO(t)

## Make figure
fig = Figure(size = (800, 700))
ax  = Axis3(fig[1, 1];
    xlabel = "GEO x (Greenwich)",
    ylabel = "GEO y (90 E)",
    zlabel = "GEO z (rotation axis)",
    title  = "SM axes in the GEO frame, $(t)",
    aspect = :data,
)

# Earth, unit radius
u = range(0, π,  length = 40)
v = range(0, 2π, length = 60)
surface!(ax,
    [sin(ui) * cos(vi) for ui in u, vi in v],
    [sin(ui) * sin(vi) for ui in u, vi in v],
    [cos(ui)           for ui in u, vi in v];
    colormap = :blues, alpha = 0.3, shading = NoShading,
)

# Draw a vector from the origin, with a label at the tip
function draw_axis!(ax, vec, color, label; length = 1.6, linewidth = 4, linestyle = :solid)
    tip = length .* vec
    lines!(ax, [0.0, tip[1]], [0.0, tip[2]], [0.0, tip[3]];
        color = color, linewidth = linewidth, linestyle = linestyle)
    scatter!(ax, [tip[1]], [tip[2]], [tip[3]]; color = color, markersize = 10)
    text!(ax, Point3f(1.1 .* tip...); text = label, color = color)
end

# GEO axes, for reference
draw_axis!(ax, [1.0, 0.0, 0.0], :gray, "GEO x"; length = 1.3, linewidth = 2, linestyle = :dash)
draw_axis!(ax, [0.0, 1.0, 0.0], :gray, "GEO y"; length = 1.3, linewidth = 2, linestyle = :dash)
draw_axis!(ax, [0.0, 0.0, 1.0], :gray, "GEO z"; length = 1.3, linewidth = 2, linestyle = :dash)

# SM axes
draw_axis!(ax, x̂, :red,   "x_SM")
draw_axis!(ax, ŷ, :green, "y_SM")
draw_axis!(ax, ẑ, :blue,  "z_SM (dipole axis)")

# Sun direction
draw_axis!(ax, s, :orange, "Sun"; linewidth = 3)

##
isdefined(Main, :server) && close(server)
server = Bonito.Server(Bonito.App(fig), "0.0.0.0", 8888)
