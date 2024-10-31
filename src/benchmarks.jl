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
              " Please run `julia $(install_script)` to do so.")
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
