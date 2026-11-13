
root  = "/mnt/data/ath099/AURORA.jl/src/WPI"
files = filter(f -> endswith(f, ".jl"), readdir(root; join=true))

## Optional: only include these functions. Set to `nothing` to include all of them.
chosen = Set([
"DentonDensity",
"FP_to_ENU",
"LLA_to_SM",
"ParticleState",
"PlasmaState",
"SM_to_ENU",
"SM_to_GEO",
"StationFieldLine",
"TsyganenkoMagneticField",
"WPI_TOF",
"WaveState",
"boris_mover",
"denton_density_model",
"dipole_field",
"dispersion_relation_whistler_branch",
"electron_density",
"group_velocity_whistler_wave",
"losscone_angle",
"magnetic_basis",
"minimal_resonance_energy",
"particle_transit",
"pitch_angle",
"point_of_resonance",
"resonance_condition",
"trace_fieldline",
"tsyganenko_field_gsm",
"velocity_from_kinetic_energy",
"wave_fields",
"wave_frequency",
"wave_launch_time",
"wave_transit",
])
#chosen = nothing

# Collect every `function name(` definition
fn_regex = r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)\("
defs = Set{String}()
for f in files, line in eachline(f)
    m = match(fn_regex, line)
    isnothing(m) || push!(defs, m.captures[1])
end
nodes = isnothing(chosen) ? defs : chosen

# An edge (a, b) means: inside the body of a, there is a call `b(`
edges = Set{Tuple{String,String}}()
for f in files
    current = nothing
    for line in eachline(f)
        m = match(fn_regex, line)
        if !isnothing(m)
            current = m.captures[1]
            continue
        end
        (isnothing(current) || current ∉ nodes) && continue
        for d in nodes
            d == current && continue
            occursin(Regex("\\b" * d * "\\("), line) && push!(edges, (current, d))
        end
    end
end

missing_names = setdiff(nodes, defs)
isempty(missing_names) || @warn "Not found as `function name(` definitions" missing_names
@show length(defs), length(edges)


outfile = "/mnt/data/ath099/AURORA.jl/callgraph.dot"
open(outfile, "w") do io
    println(io, "digraph G { rankdir=TB; node [shape=box];")
    for name in sort(collect(nodes))          # declare every node, so isolated ones also show
        println(io, "  \"$name\";")
    end
    for (a, b) in sort(collect(edges))
        println(io, "  \"$a\" -> \"$b\";")
    end
    println(io, "}")
end

#===
Render  in the terminal:

dot -Tpng /mnt/data/ath099/AURORA.jl/callgraph.dot -o /mnt/data/ath099/AURORA.jl/callgraph.png
===#
