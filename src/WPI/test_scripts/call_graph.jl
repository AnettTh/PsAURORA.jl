##
files = filter(f -> endswith(f, ".jl"), readdir("src/WPI"; join=true))
# 1. Collect function names (only `function name(` definitions)
defs = Set{String}()
for f in files, line in eachline(f)
    m = match(r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)\(", line)
    isnothing(m) || push!(defs, m.captures[1])
end
# 2. For each function body, find calls to the other collected names
edges = Set{Tuple{String,String}}()
for f in files
    current = nothing
    for line in eachline(f)
        m = match(r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)\(", line)
        if !isnothing(m)
            current = m.captures[1]
            continue
        end
        isnothing(current) && continue
        for d in defs
            d == current && continue
            occursin(Regex("\\b" * d * "\\("), line) && push!(edges, (current, d))
        end
    end
end
# 3. Write a Graphviz file
open("callgraph.dot", "w") do io
    println(io, "digraph G { rankdir=LR; node [shape=box];")
    for (a, b) in sort(collect(edges))
        println(io, "  \"$a\" -> \"$b\";")
    end
    println(io, "}")
end


#===

Call-graph with selected functions only

===#

## All functions to choose from
files = filter(f -> endswith(f, ".jl"), readdir("/mnt/data/ath099/AURORA.jl/src/WPI"; join=true))

defs = Set{String}()
for f in files, line in eachline(f)
    m = match(r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)\(", line)
    isnothing(m) || push!(defs, m.captures[1])
end

foreach(println, sort(collect(defs)))
##
small_set = Set([
    "DentonDensity",
    "ParticleState",
    "PlasmaState",
    "TsyganenkoMagneticField",
    "WPI_TOF",
    "WaveState",
    "boris_mover",
    "denton_density_model",
    "DipoleMagneticField",
    "dipole_field",
    "dipole_tilt",
    "dispersion_relation_whistler_branch",
    "electron_density",
    "find_L_shell",
    "find_R_max",
    "find_footpoint",
    "group_velocity_whistler_wave",
    "losscone_angle",
    "magnetic_basis",
    "particle_transit",
    "point_on_line",
    "resonance_condition",
    "resonance_latitude",
    "trace_fieldline",
    "tsyganenko_field_gsm",
    "velocity_from_kinetic_energy",
    "wave_fields",
    "wave_frequency",
    "wave_launch_time",
    "wave_transit",
])

##
chosen = Set([
    "Cartesian",
    "DentonDensity",
    "Dipolar",
    "FP_to_ENU",
    "LLA_to_SM",
    "ParticleState",
    "PlasmaState",
    "SM_to_ENU",
    "SM_to_GEO",
    "TsyganenkoMagneticField",
    "WPI_TOF",
    "WaveState",
    "axis_SM",
    "boris_mover",
    "denton_density_model",
    "dipole_axis_GEO",
    "dipole_field",
    "dipole_tilt",
    "dispersion_relation_whistler_branch",
    "electron_density",
    "find_L_shell",
    "find_R_max",
    "find_footpoint",
    "group_velocity_whistler_wave",
    "longitude_to_MLT",
    "losscone_angle",
    "magnetic_basis",
    "particle_transit",
    "point_on_line",
    "resonance_condition",
    "resonance_latitude",
    "solar_coefficients",
    "subsolar_point",
    "sun_unit_vector_GEO",
    "trace_fieldline",
    "tsyganenko_field_gsm",
    "velocity_from_kinetic_energy",
    "wave_fields",
    "wave_frequency",
    "wave_launch_time",
    "wave_transit",
])

edges = Set{Tuple{String,String}}()
for f in files
    current = nothing
    for line in eachline(f)
        m = match(r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)\(", line)
        if !isnothing(m)
            current = m.captures[1]
            continue
        end
        (isnothing(current) || current ∉ chosen) && continue     # caller must be chosen
        for d in chosen                                          # callee must be chosen
            d == current && continue
            occursin(Regex("\\b" * d * "\\("), line) && push!(edges, (current, d))
        end
    end
end

missing_names = setdiff(chosen, defs)
isempty(missing_names) || @warn "Not found as `function name(` definitions" missing_names

outfile = "/mnt/data/ath099/AURORA.jl/callgraph.dot"
open(outfile, "w") do io
    println(io, "digraph G { rankdir=TB; node [shape=box];")
    for name in sort(collect(chosen))          # declare every node, so isolated ones also show
        println(io, "  \"$name\";")
    end
    for (a, b) in sort(collect(edges))
        println(io, "  \"$a\" -> \"$b\";")
    end
    println(io, "}")
end
@show length(edges)
