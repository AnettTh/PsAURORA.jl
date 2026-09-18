using AURORA

# NOTE: For several energies
## Define the particles
μ = -0.99       # Almost field-aligned
L = 6.0
r0 = [L*RE, 0.0, 0.0]

E_grid = range(1e3, 10000e3, length=100)
particles = [ParticleState(E, μ, r0, dipole_field) for E in E_grid]

## Define the plasma
λ_grid = range(0.0, deg2rad(50), length=500)
n_e0 = 1.8e7        # 18/cc, from Hsieh 2022

plasma = PlasmaParameters(λ_grid, n_e0, dipole_field, L)


## Define the wave
ωs = range(plasma.Ω_e[1]*0.2, plasma.Ω_e[1]*0.4, length=5)


## Get parallel velocity
v_parallel_grid = [parallel_velocity(p, λ_grid) for p in particles]  # [n_E][n_λ]


## Get group velocity
v_g_set = [group_velocity_whistler_wave(ω; Ω_e=plasma.Ω_e[1], ω_pe=plasma.ω_pe[1]) for ω in ωs]


## Calculate TOF
TOF_matrix = zeros(length(E_grid), length(ωs))  # [n_E × n_ω]
TOF_simple = zeros(length(E_grid), length(ωs))

for (i, p) in enumerate(particles)
    TOF_matrix[i, :] = WPI_TOF_field_dependent(ωs, plasma, p, t0)
    v_parallel_lc = p.v * cos(p.α_lc)
    for (j, v_g) in enumerate(v_g_set)
        TOF_simple[i, j] = WPI_TOF(deg2rad(20), v_g, v_parallel_lc, p.L*RE, 0, p.E_eV, p.α_lc)
    end
end


## Make Figure
fig = Figure()
ax = Axis(fig[1,1], xlabel="time-of-flight [s]", ylabel="E [keV]", yscale=log10)

ax.yticks = [1, 10, 100, 1000]
ax.ytickformat = values -> ["$(Int(v))" for v in values]

ω_norm = (ωs .- minimum(ωs)) ./ (maximum(ωs) - minimum(ωs))


colors = [:blue, :red, :green, :orange, :purple]  # one per ω

for (i, ω) in enumerate(ωs)
    lines!(ax, TOF_matrix[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:solid,
        label="ω = $(round(ω/1e3, digits=1)) kHz"
    )
    lines!(ax, TOF_simple[:, i], E_grid ./ 1e3;
        color=colors[i],
        linestyle=:dash
    )
end

axislegend(ax, position=:rt)

#xlims!(ax, 0.6, 0.65)
#ylims!(ax, 100, 1000)

fig
