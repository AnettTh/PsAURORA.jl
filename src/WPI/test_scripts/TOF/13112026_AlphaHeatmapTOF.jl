using AURORA
using CairoMakie
using Profile
using ProgressMeter
using JLD2


## Define magnetic field

## Define the particles

## Define the plasma

# Define the wave

## Calculate TOF
TOF = zeros(length(μ0_grid), length(ωs))
progress = Progress(length(particle))

Threads.@threads :greedy for i in eachindex(particle)
    TOF[i, :] = WPI_TOF(
        ωs,
        particle[i],
        plasma;
        wave=wave,
    )
    next!(progress)
end


## Make Figure
fig = Figure(size=(700,500))
ax = Axis(
    fig[1,1],
    xlabel="ω/Ωₑ",
    ylabel="α [°]",
    title="Time-of-flight using "
    * if use_boris "boris mover \n" else "integration \n " end
    * "E=$(E / 1e3) keV"
    * ", MLT $(round(longitude_to_MLT(ϕ_eq; degrees=false), digits = 1))"
    * ", L=$(L)",
)

hm = heatmap!(ax, ωs ./ Ω_e_eq, rad2deg.(acos.(-μ0_grid)), TOF'; colormap=:turbo)
Colorbar(fig[1,2], hm; label="TOF [s]")
##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

save("src/WPI/test_scripts/TOF/figures/02112026_TOF_alpha_heatmap.png", fig)
