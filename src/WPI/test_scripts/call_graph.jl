##
files = filter(f -> endswith(f, ".jl"), readdir("src/WPI"; join=true))

## 1. Collect function names (only `function name(` definitions)
defs = Set{String}()
for f in files, line in eachline(f)
    m = match(r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)\(", line)
    isnothing(m) || push!(defs, m.captures[1])
end

## 2. For each function body, find calls to the other collected names
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

## 3. Write a Graphviz file
open("callgraph.dot", "w") do io
    println(io, "digraph G { rankdir=LR; node [shape=box];")
    for (a, b) in sort(collect(edges))
        println(io, "  \"$a\" -> \"$b\";")
    end
    println(io, "}")
end
