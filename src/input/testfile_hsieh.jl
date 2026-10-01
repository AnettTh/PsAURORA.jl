using AURORA

##
L = 4.5
n_e = 18e6
n_h = 0.005 * n_e

##
V_t_parallel = 0.05 * c₀
V_t_perp = 0.1 * c₀

##
ρ = 0.2
β = 0.3

##
K_grid = range(1e3 * eV_in_J, 300e3 * eV_in_J, length=300)  # [J]
α_grid = range(1e-3, π/2, length=300)                          # [rad], avoid 0

K_2d = K_grid'              # 1 × n_K
α_2d = collect(α_grid)      # n_α × 1

F_Kα = [subtracted_bimaxwellian_Hsieh(K, α, n_h, V_t_parallel, V_t_perp, ρ, β)
        for α in α_grid, K in K_grid]


##
F_Kα_log = @. ifelse(F_Kα ≤ 0, NaN, log10(F_Kα))
F_max    = maximum(filter(isfinite, vec(F_Kα_log)))
clims    = (F_max - 6, F_max)

fig = Figure(size=(700, 500))
ax  = Axis(fig[1, 1];
    xlabel = "K [keV]",
    ylabel = "α [°]",
    title  = "Hsieh 2022 F_EQ(K, α)"
)

hm = heatmap!(ax,
    collect(K_grid) ./ eV_in_J ./ 1e3,
    rad2deg.(α_grid),
    F_Kα_log;
    colormap   = :inferno,
    colorrange = clims
)

Colorbar(fig[1, 2], hm2; label="log₁₀(F_EQ)", tellheight=false)
ylims!(0, 15)
fig
