using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Scratch: get_scratch!
using UUIDs: UUID
using NaturalSort: natural
using ZipArchives: ZipReader, zip_names, zip_readentry
import Git

function install_baharev_benchmarks()
    # UUID for FeedbackArcSets.
    uuid = UUID("6c3ede71-d29b-41ca-966d-1d2ca331f31c")
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
    open(index_file, "a") do io
        foreach(sort(index, lt = natural)) do entry
            println(io, entry)
        end
    end
end

install_baharev_benchmarks()
