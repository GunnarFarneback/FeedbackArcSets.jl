using Scratch: @get_scratch!
using Graphs: SimpleDiGraph, add_edge!

struct BaharevBenchmark
    name::String
    graph::SimpleDiGraph
    feedback_arc_set::Vector{Tuple{Int, Int}}
    cycles::Vector{Vector{Int}}
end

function Base.show(io::IO, bench::BaharevBenchmark)
    print(io, "BaharevBenchmark: $(bench.name)")
end

"""
    baharev_benchmark()

Return all 4468 benchmarks used in "An Exact Method for the Minimum
Feedback Arc Set Problem" by Baharev et al. The return value is a
vector of structs with the fields:
* `name`: The name of the benchmark case.
* `graph`: The graph in the form of a `Graphs.SimpleDiGraph`.
* `feedback_arc_set`: One minimal feedback arc set represented as a
  vector of tuples.
* `cycles`: A set of cycles sufficient to prove the lower bound of the
  feedback arc set size when used as constraints in an Integer Linear
  Program.

Notes:
* Compared to the original sources, all vertex numbers have been
  increased by one to match Julia conventions.
* The benchmark case `"erdos_renyi_n_85_c_5_seed_24"` provides a
  non-minimal feedback arc set.


    baharev_benchmark(index::Integer)

Return a single benchmark with the given `index` in the full set of
benchmarks.

    baharev_benchmark(indices::AbstractVector{<:Integer})

Return a selection of benchmarks with the given `indices`.

    baharev_benchmark(name::AbstractString)

Return a single benchmark with the given `name`.

    baharev_benchmark(names::AbstractVector{<:AbstractString})

Return a selection of the benchmarks with the given `names`.

    baharev_benchmark(predicate::Function)

Return a selection of benchmarks for which `predicate(name)` returns
true.

# Examples:
    baharev_benchmark(13)
    baharev_benchmark(1:5)
    baharev_benchmark("Imase_Itoh_n_100_d_3")
    baharev_benchmark(["tournament_n_28_seed_1", "erdos_renyi_n_35_c_5_seed_1"])
    baharev_benchmark(startswith("de_Bruijn"))
"""
function baharev_benchmark(selection::Union{AbstractVector{<:AbstractString},
                                            AbstractVector{<:Integer},
                                            Function})
    data_dir = @get_scratch!("baharev")
    index_file = joinpath(data_dir, "index.txt")

    if !isfile(index_file)
        install_script = abspath(@__DIR__, "..", "benchmarks",
                                 "install_benchmark_data.jl")
        error("The Baharev benchmark data (79 MB) has not been installed." *
              " Please run `julia $(install_script)` to install it.")
    end

    if selection isa Function
        names = filter(selection, readlines(index_file))
    elseif eltype(selection) <: Integer
        names = readlines(index_file)[selection]
    else
        names = selection
    end

    benchmarks = BaharevBenchmark[]

    for name in names
        edges = _read_baharev_file(data_dir, name, ".edges")
        mfes = _read_baharev_file(data_dir, name, ".mfes")
        cycles = _read_baharev_file(data_dir, name, ".cycles")
        num_vertices = maximum(reduce(vcat, edges))
        graph = SimpleDiGraph(num_vertices)
        for edge in edges
            add_edge!(graph, edge...)
        end
        push!(benchmarks,
              BaharevBenchmark(name, graph, Tuple.(mfes), cycles))
    end

    return benchmarks
end

baharev_benchmark(n::Integer) = only(baharev_benchmark([n]))
baharev_benchmark(name::AbstractString) = only(baharev_benchmark([name]))
baharev_benchmark() = baharev_benchmark(Returns(true))

function _read_baharev_file(data_dir, name, suffix)
    return [parse.(Int, split(line)) .+ 1
            for line in eachline(joinpath(data_dir, name * suffix))]
end

struct DasdanBenchmark
    name::String
    graph::SimpleDiGraph
    directory::String
end

function Base.show(io::IO, bench::DasdanBenchmark)
    print(io, "DasdanBenchmark: $(bench.name)")
end

"""
    dasdan_benchmark()

Return all 141 benchmark graphs from the repository
https://github.com/alidasdan/graph-benchmarks. The return value is a
vector of structs with the fields:
* `directory`: The name of the subdirectory containing the graph.
* `name`: The name of the graph.
* `graph`: The graph in the form of a `Graphs.SimpleDiGraph`.

Notes:

* The original sources also contain an "arc weight" and "transit time"
  for each edge in the graphs. Those are not considered by this
  function.


    dasdan_benchmark(index::Integer)

Return a single benchmark with the given `index` in the full set of
benchmarks.

    dasdan_benchmark(indices::AbstractVector{<:Integer})

Return a selection of benchmarks with the given `indices`.

    dasdan_benchmark(name::AbstractString)

Return a single benchmark with the given `name`.

    dasdan_benchmark(names::AbstractVector{<:AbstractString})

Return a selection of the benchmarks with the given `names`.

    dasdan_benchmark(predicate::Function)

Return a selection of benchmarks for which `predicate(name)` returns
true.

# Examples:
    dasdan_benchmark(13)
    dasdan_benchmark(1:5)
    dasdan_benchmark("dsip")
    dasdan_benchmark(["sample", "r05"])
    dasdan_benchmark(startswith("ibm"))
"""
function dasdan_benchmark(selection::Union{AbstractVector{<:AbstractString},
                                           AbstractVector{<:Integer},
                                           Function})
    data_dir = @get_scratch!("dasdan")
    index_file = joinpath(data_dir, "index.txt")

    if !isfile(index_file)
        install_script = abspath(@__DIR__, "..", "benchmarks",
                                 "install_benchmark_data.jl")
        error("The Dasdan benchmark data (197 MB) has not been installed." *
              " Please run `julia $(install_script)` to install it.")
    end

    name_to_dir = Dict(last(split(line)) => first(split(line))
                       for line in eachline(index_file))
    all_names = [last(split(line)) for line in eachline(index_file)]
    if selection isa Function
        names = filter(selection, all_names)
    elseif eltype(selection) <: Integer
        names = all_names[selection]
    else
        names = selection
    end

    benchmarks = DasdanBenchmark[]

    for name in names
        local graph
        for line in eachline(joinpath(data_dir, "$(name).d"))
            isempty(line) && continue
            parts = split(line)
            if first(parts) == "p"
                graph = SimpleDiGraph(parse(Int, parts[3]))
            elseif first(parts) == "a"
                add_edge!(graph, parse.(Int, parts[[2, 3]])...)
            end
        end
        push!(benchmarks,
              DasdanBenchmark(name, graph, name_to_dir[name]))
    end

    return benchmarks
end

dasdan_benchmark(n::Integer) = only(dasdan_benchmark([n]))
dasdan_benchmark(name::AbstractString) = only(dasdan_benchmark([name]))
dasdan_benchmark() = dasdan_benchmark(Returns(true))

struct SnapBenchmark
    name::String
    graph::SimpleDiGraph
end

function Base.show(io::IO, bench::SnapBenchmark)
    print(io, "SnapBenchmark: $(bench.name)")
end

"""
    snap_benchmark()

Return 16 directed graphs from the [Stanford Large Network Dataset
Collection](http://snap.stanford.edu/data/). This selection is used as
benchmark in "Improvement and Evaluation of a Heuristic Method for the
Minimal Feedback Arc Set Problem" by Pustoslemšek et al.

The return value is a vector of structs with the fields:
* `name`: The name of the graph.
* `graph`: The graph in the form of a `Graphs.SimpleDiGraph`.

Notes:

* The vertex numbers from the original sources are used unchanged
  unless they start with 0, in which case all vertex numbers are
  increased by one. In some cases this means that a few vertices are
  isolated without any edges.

    snap_benchmark(index::Integer)

Return a single benchmark with the given `index` in the full set of
benchmarks.

    snap_benchmark(indices::AbstractVector{<:Integer})

Return a selection of benchmarks with the given `indices`.

    snap_benchmark(name::AbstractString)

Return a single benchmark with the given `name`.

    snap_benchmark(names::AbstractVector{<:AbstractString})

Return a selection of the benchmarks with the given `names`.

    snap_benchmark(predicate::Function)

Return a selection of benchmarks for which `predicate(name)` returns
true.

# Examples:
    snap_benchmark(13)
    snap_benchmark(1:5)
    snap_benchmark("soc-Epinions1.txt")
    snap_benchmark(["email-EuAll.txt", "web-Stanford.txt"])
    snap_benchmark(startswith("p2p-Gnutella"))
"""
function snap_benchmark(selection::Union{AbstractVector{<:AbstractString},
                                         AbstractVector{<:Integer},
                                         Function})
    data_dir = @get_scratch!("snap")
    index_file = joinpath(data_dir, "index.txt")

    if !isfile(index_file)
        install_script = abspath(@__DIR__, "..", "benchmarks",
                                 "install_benchmark_data.jl")
        error("The Snap benchmark data (92 MB) has not been installed." *
              " Please run `julia $(install_script)` to install it.")
    end

    all_names = readlines(index_file)
    if selection isa Function
        names = filter(selection, all_names)
    elseif eltype(selection) <: Integer
        names = all_names[selection]
    else
        names = selection
    end

    benchmarks = SnapBenchmark[]

    for name in names
        edges = Vector{Int}[]
        for line in eachline(joinpath(data_dir, "$(name).txt"))
            startswith(line, "#") && continue
            push!(edges, parse.(Int, split(line)))
        end
        first_vertex, last_vertex = extrema(reduce(vcat, edges))
        @assert first_vertex >= 0
        offset = (first_vertex == 0) ? 1 : 0
        graph = SimpleDiGraph(last_vertex + offset)
        for (src, dst) in edges
            add_edge!(graph, src + offset, dst + offset)
        end
        push!(benchmarks, SnapBenchmark(name, graph))
    end

    return benchmarks
end

snap_benchmark(n::Integer) = only(snap_benchmark([n]))
snap_benchmark(name::AbstractString) = only(snap_benchmark([name]))
snap_benchmark() = snap_benchmark(Returns(true))

struct SmallGraphsBenchmark
    graph::SimpleDiGraph
end

function Base.show(io::IO, bench::SmallGraphsBenchmark)
    v = nv(bench.graph)
    e = ne(bench.graph)
    print(io, "SmallGraphsBenchmark: $v vertices, $e edges")
end

"""
    small_graphs_benchmark()

Return all 569, up to isomorphism, strongly connected graphs with 2 to
8 edges.

The return value is a vector of structs with the field:
* `graph`: The graph in the form of a `Graphs.SimpleDiGraph`.

    small_graphs_benchmark(index::Integer)

Return a single benchmark with the given `index` in the full set of
benchmarks.

    small_graphs_benchmark(indices::AbstractVector{<:Integer})

Return a selection of benchmarks with the given `indices`.

# Examples:
    small_graphs_benchmark(13)
    small_graphs_benchmark(1:5)
"""
function small_graphs_benchmark(selection::AbstractVector{<:Integer})
    data_dir = @get_scratch!("small_graphs")
    index_file = joinpath(data_dir, "index.txt")
    data_file = joinpath(data_dir, "graphs.bin")

    if !isfile(index_file)
        install_script = abspath(@__DIR__, "..", "benchmarks",
                                 "install_benchmark_data.jl")
        error("The small_graphs benchmark data (5 kB) has not been installed." *
              " Please run `julia $(install_script)` to install it." *
              " It may take some time to generate it.")
    end

    encoded_graphs = reinterpret(UInt64, read(data_file))
    return decode_small_graph.(encoded_graphs[selection])
end

small_graphs_benchmark(n::Integer) = only(small_graphs_benchmark([n]))
small_graphs_benchmark() = small_graphs_benchmark(1:569)

function decode_small_graph(x)
    num_vertices = ceil(Int, sqrt(log2(x)))
    adj = reshape((x .>> (0:(num_vertices^2 -1))) .& 1,
                  num_vertices, num_vertices)
    return SimpleDiGraph(adj)
end
