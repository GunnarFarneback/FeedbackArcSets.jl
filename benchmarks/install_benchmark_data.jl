using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Scratch: get_scratch!
using UUIDs: UUID
using NaturalSort: natural
using ZipArchives: ZipReader, zip_names, zip_readentry
using Inflate: inflate_gzip
using Downloads: download
import Git
using REPL.TerminalMenus: MultiSelectMenu, request
using NautyGraphs: NautyDiGraph, ghash, canonize!
using Graphs: SimpleDiGraph, adjacency_matrix, is_strongly_connected,
              cycle_digraph
using Combinatorics: combinations

# UUID for FeedbackArcSets.
const uuid = UUID("6c3ede71-d29b-41ca-966d-1d2ca331f31c")

function install_baharev_benchmarks()
    data_dir = get_scratch!(uuid, "baharev")
    index_file = joinpath(data_dir, "index.txt")
    if isfile(index_file)
        @info "The Baharev benchmarks are already installed."
        return
    end
    index = String[]
    mktempdir() do tmpdir
        url = "https://github.com/baharev/sdopt-tearing.git"
        run(`$(Git.git()) -C $tmpdir clone $url`)
        zipdir = joinpath(tmpdir, "sdopt-tearing", "data",
                          "benchmark_mfes")
        zipfiles = filter(endswith(".zip"), readdir(zipdir, join = true))
        for zipfile in zipfiles
            archive = ZipReader(read(zipfile))
            for entry in zip_names(archive)
                if last(splitext(entry)) in (".edges", ".mfes", ".cycles")
                    filename = basename(entry)
                    write(joinpath(data_dir, filename),
                          zip_readentry(archive, entry))
                    if endswith(filename, ".edges")
                        push!(index, first(splitext(filename)))
                    end
                end
            end
        end
    end
    open(index_file, "w") do io
        foreach(sort(index, lt = natural)) do entry
            println(io, entry)
        end
    end
end

function install_dasdan_benchmarks()
    data_dir = get_scratch!(uuid, "dasdan")
    index_file = joinpath(data_dir, "index.txt")
    if isfile(index_file)
        @info "The Dasdan benchmarks are already installed."
        return
    end
    index = String[]
    mktempdir() do tmpdir
        url = "https://github.com/alidasdan/graph-benchmarks.git"
        for (root, dirs, files) in walkdir(tmpdir)
            graph_files = filter(name -> (endswith(name, ".d")
                                          || endswith(name, ".d.gz")),
                                 files)
            if !isempty(graph_files)
                for file in graph_files
                    source_file = joinpath(root, file)
                    target_file = joinpath(data_dir,
                                           replace(file, ".d.gz" => ".d"))
                    name = first(split(file, "."))
                    dir = basename(root)
                    if endswith(file, ".gz")
                        write(target_file, inflate_gzip(source_file))
                    else
                        write(target_file, read(source_file))
                    end
                    push!(index, "$dir $name")
                end
            end
        end
    end
    open(index_file, "w") do io
        foreach(sort(index, lt = natural)) do entry
            println(io, entry)
        end
    end
end

function install_snap_benchmarks()
    data_dir = get_scratch!(uuid, "snap")
    index_file = joinpath(data_dir, "index.txt")
    if isfile(index_file)
        @info "The SNAP benchmarks are already installed."
        return
    end
    index = String[]
    mktempdir() do tmpdir
        urls = ["http://snap.stanford.edu/data/wiki-Vote.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella08.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella09.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella06.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella05.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella04.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella25.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella24.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella30.txt.gz",
                "http://snap.stanford.edu/data/p2p-Gnutella31.txt.gz",
                "http://snap.stanford.edu/data/soc-Epinions1.txt.gz",
                "http://snap.stanford.edu/data/email-EuAll.txt.gz",
                "http://snap.stanford.edu/data/web-NotreDame.txt.gz",
                "http://snap.stanford.edu/data/web-Stanford.txt.gz",
                "http://snap.stanford.edu/data/soc-Slashdot0811.txt.gz",
                "http://snap.stanford.edu/data/soc-Slashdot0902.txt.gz"]

        for url in urls
            filename_gz = joinpath(tmpdir, basename(url))
            filename, _ = splitext(basename(url))
            name, _ = splitext(filename)
            target_file = joinpath(data_dir, filename)
            download(url, filename_gz)
            write(target_file, inflate_gzip(filename_gz))
            push!(index, name)
        end
    end
    open(index_file, "w") do io
        foreach(sort(index, lt = natural)) do entry
            println(io, entry)
        end
    end
end

# Adjacency matrix for a single edge from v to w.
function edge_matrix(num_vertices, v, w)
    a = zeros(Int, num_vertices, num_vertices)
    a[v, w] = 1
    return a
end

# Generate small strongly connected graphs with given number of
# vertices and edges.
function small_graphs(num_vertices, num_edges)
    edges = [edge_matrix(num_vertices, v, w)
             for v in 1:num_vertices, w in 1:num_vertices
             if v != w]
    # We can without loss of generality assume that the edge from
    # vertex 1 to 2 is part of the graph. This is not necessary to do
    # but saves lot of time for the higher numbers of vertices and
    # edges.
    fixed_edges = popfirst!(edges)
    num_fixed_edges = 1
    # Furthermore any strongly connected graph with more than two
    # edges must contain at least one combination of edges a -> b and
    # b -> c. Thus we can without loss of generality also assume that
    # 2 -> 3 is part of the graph.
    if num_vertices > 2 && num_edges > 2
        fixed_edges += popat!(edges, num_vertices)
        num_fixed_edges += 1
    end
    graphs = Dict{UInt64, NautyDiGraph}()
    for combination in combinations(edges, num_edges - num_fixed_edges)
        adjacency_matrix = fixed_edges + sum(combination)
        g = NautyDiGraph(adjacency_matrix)
        h = ghash(g)
        if !haskey(graphs, h)
            canonize!(g)
            graphs[h] = g
        end
    end
    strongly_connected_graphs = SimpleDiGraph[]
    for graph in values(graphs)
        g = SimpleDiGraph(adjacency_matrix(graph))
        if is_strongly_connected(g)
            push!(strongly_connected_graphs, g)
        end
    end
    return strongly_connected_graphs
end

function encode_small_graph(graph)
    a = vec(adjacency_matrix(graph))
    return reduce(|, UInt64.(a) .<< (0:(length(a) - 1)))
end

function install_small_graphs_benchmarks()
    data_dir = get_scratch!(uuid, "small_graphs")
    index_file = joinpath(data_dir, "index.txt")
    if isfile(index_file)
        @info "The small_graphs benchmarks are already installed."
        return
    end
    data_file = joinpath(data_dir, "graphs.bin")
    graph_counts = Dict(num_edges => 0 for num_edges in 2:8)
    open(data_file, "w") do io
        for num_edges in 2:8
            for num_vertices in 2:num_edges
                if (num_vertices, num_edges) == (8, 8)
                    # It takes unnecessarily much time to search for this
                    # single graph.
                    graphs = [cycle_digraph(8)]
                else
                    graphs = small_graphs(num_vertices, num_edges)
                end
                graph_counts[num_edges] += length(graphs)
                for graph in graphs
                    write(io, encode_small_graph(graph))
                end
            end
        end
    end

    write(index_file,
          """
          Info: The small graphs are all stored in binary form in graphs.bin,
                one UInt64 value each, where the bits form the adjacency
                matrices.
          """)

    # Sanity check that we have the correct number of graphs.
    # Cf. https://oeis.org/A350752.
    if [graph_counts[i] for i in 2:8] != [1, 1, 3, 6, 25, 91, 442]
        error("Internal error. An incorrect number of small graphs were generated.")
    end
end

function install_benchmarks()
    benchmarks = [["baharev", 79, install_baharev_benchmarks],
                  ["dasdan", 197, install_dasdan_benchmarks],
                  ["snap", 92, install_snap_benchmarks],
                  ["small_graphs", 0.005, install_small_graphs_benchmarks]]

    installed_benchmarks = String[]
    uninstalled_benchmarks = []
    for benchmark in benchmarks
        name = benchmark[1]
        if isfile(joinpath(get_scratch!(uuid, name), "index.txt"))
            push!(installed_benchmarks, name)
        else
            push!(uninstalled_benchmarks, benchmark)
        end
    end

    if !isempty(installed_benchmarks)
        println("The following benchmark sets are installed:")
        for name in installed_benchmarks
            println(name)
        end
        println()
    end

    if !isempty(uninstalled_benchmarks)
        println("Choose which benchmark sets to install.")
        options = [string(benchmark[1], " (", benchmark[2], " MB)")
                   for benchmark in uninstalled_benchmarks]
        selected = request(MultiSelectMenu(options))
        println()
        for i in sort(collect(selected))
            name, _, installer = uninstalled_benchmarks[i]
            println("Installing $(name) benchmark data.")
            installer()
        end
        isempty(selected) || println()
    end

    println("To remove benchmark data, manually delete the corresponding directory in\n",
            dirname(get_scratch!(uuid, "baharev")))
end

install_benchmarks()
