using Test: @testset, @test, @test_throws
using FeedbackArcSets: find_feedback_arc_set, dfs_feedback_arc_set,
                       greedy_feedback_arc_set, pagerank_feedback_arc_set,
                       is_feedback_arc_set
using Graphs: path_digraph, cycle_digraph, complete_digraph, SimpleDiGraph,
              add_edge!, rem_edge!
using GoGameGraphs: go_game_graph
import HiGHS

heuristic_methods = [dfs_feedback_arc_set,
                     greedy_feedback_arc_set,
                     g -> greedy_feedback_arc_set(g, randomize = false),
                     pagerank_feedback_arc_set,
                     g -> pagerank_feedback_arc_set(g, num_iterations = 1)]

# Paths are acyclic and all methods should return an empty feedback arc set.
@testset "path graphs" begin
    for n = 0:10
        g = path_digraph(n)
        exact = find_feedback_arc_set(g; log_level = 0, split = false)
        @test exact.lower_bound == length(exact.feedback_arc_set) == 0
        for f in heuristic_methods
            @test isempty(f(g))
        end
    end
end

# Cycle graphs have a single Hamiltonian cycle and one arbitrary edge
# is a minimal feedback arc set. All methods should find that.
@testset "cycle graphs" begin
    for n = 2:10
        g = cycle_digraph(n)
        exact = find_feedback_arc_set(g; log_level = 0, split = false)
        @test exact.lower_bound == length(exact.feedback_arc_set) == 1
        for f in heuristic_methods
            @test length(f(g)) == 1
        end
    end
end

# Complete digraphs have length two cycles between each pair of
# vertices. Obviously one edge from each of these must be part of the
# minimal feedback arc set and if they are chosen with a consistent
# direction, an acyclic lattice remains. It is not immediately obvious
# that all heuristic methods should find this but apparently they do,
# so test for it.
@testset "complete graphs" begin
    for n = 2:10
        g = complete_digraph(n)
        exact = find_feedback_arc_set(g, log_level = 0, split = false)
        N = n * (n - 1) ÷ 2
        @test exact.lower_bound == length(exact.feedback_arc_set) == N
        @test is_feedback_arc_set(g, exact.feedback_arc_set)
        for f in heuristic_methods
            arc_set = f(g)
            @test length(arc_set) == N
            @test is_feedback_arc_set(g, arc_set)
        end
    end
end

# Tournament graphs have one edge between each pair of vertices. The
# construction below has no particular properties. It's only meant to
# produce a deterministic and not entirely trivial tournament graph.
function tournament_graph(n)
    g = SimpleDiGraph(n)
    for i = 1:n
        for j = (i + 1):n
            if mod(i^2 + j^3, 5) > 2
                add_edge!(g, i, j)
            else
                add_edge!(g, j, i)
            end
        end
    end
    return g
end

# Larger tournament graphs are demanding to find minimal feedback arc
# sets in but for the tested sizes it can be done quickly, at least
# with the construction used here.
#
# The heuristic methods are only tested to verify that they produce
# valid feedback arc sets, even though they may find optimal results
# for the smallest sizes.
@testset "tournament graphs" begin
    for n = 3:24
        g = tournament_graph(n)
        exact = find_feedback_arc_set(g, log_level = 0, split = false)
        @test exact.lower_bound == length(exact.feedback_arc_set)
        @test is_feedback_arc_set(g, exact.feedback_arc_set)
        for f in heuristic_methods
            arc_set = f(g)
            @test length(arc_set) >= exact.lower_bound
            @test is_feedback_arc_set(g, arc_set)
        end
    end
end

# Go game graphs are sparse, have lots of cycles, and a small
# radius. Empirically they are relatively easy to find minimal
# feedback arc sets in, so testing can be done on more than tiny
# instances. The largest graph in these tests has 981 vertices and
# 5258 edges, yet the exact solution can be found within a few
# seconds.
#
# The heuristic methods are only tested to verify that they return
# valid feedback arc sets, even though they may find optimal results
# for the smallest instances.
@testset "go game graphs" begin
    for id in [1, 3, 7, 11, 12, 13, 15, 95, 127, 1111, 33983]
        g = go_game_graph(id)
        exact = find_feedback_arc_set(g, log_level = 0, split = false)
        @test exact.lower_bound == length(exact.feedback_arc_set)
        @test is_feedback_arc_set(g, exact.feedback_arc_set)
        for f in heuristic_methods
            arc_set = f(g)
            @test length(arc_set) >= exact.lower_bound
            @test is_feedback_arc_set(g, arc_set)
        end
    end
end

# Do some tests on residual graphs after removing the feedback arc
# set. None of the algorithms should find a non-empty feedback arc set
# (even though that would be valid).
@testset "acyclic graphs" begin
    for id in [1, 3, 7, 11, 12, 13, 15, 95, 127, 1111, 33983]
        g = go_game_graph(id)
        for (i, j) in find_feedback_arc_set(g, log_level = 0, split = false).feedback_arc_set
            rem_edge!(g, i, j)
        end
        for split in (false, true)
            exact = find_feedback_arc_set(g; log_level = 0, split)
            @test exact.lower_bound == 0
            @test isempty(exact.feedback_arc_set)
        end
        for f in heuristic_methods
            @test isempty(f(g))
        end
    end
end

# Compare supported solvers. Test the JuMP interface with HiGHS.
@testset "solvers" begin
    graphs = vcat(go_game_graph.([1, 3, 7, 11, 15, 95, 127, 1111, 33983]),
                  tournament_graph.(3:24))
    solvers = [("cbc", Dict()),
               ("highs", Dict()),
               ("jump", Dict("optimizer" => HiGHS.Optimizer))]
    for g in graphs
        reference = -1
        for (solver, solver_options) in solvers
            result = find_feedback_arc_set(g; log_level = 0, split = false,
                                           solver, solver_options)
            if reference < 0
                reference = result.lower_bound
            else
                @test result.lower_bound == reference
            end
            @test is_feedback_arc_set(g, result.feedback_arc_set)
            @test length(result.feedback_arc_set) == reference
        end
    end
    @test_throws ErrorException find_feedback_arc_set(go_game_graph(1),
                                                      solver = "",
                                                      split = false)
end

@testset "self-loops" begin
    g = SimpleDiGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 1)
    add_edge!(g, 2, 1)
    add_edge!(g, 1, 1)
    add_edge!(g, 3, 3)
    @test_throws ErrorException find_feedback_arc_set(g, split = false)
    @test find_feedback_arc_set(g, self_loops = "ignore", log_level = 0,
                                split = false).feedback_arc_set == [(1, 2)]
    @test find_feedback_arc_set(g, self_loops = "ignore", log_level = 0,
                                split = false).lower_bound == 1
    @test Set(find_feedback_arc_set(g, self_loops = "include", log_level = 0,
                                    split = false).feedback_arc_set) == Set([(1, 2), (1, 1), (3, 3)])
    @test find_feedback_arc_set(g, self_loops = "include", log_level = 0,
                                split = false).lower_bound == 3
    @test_throws ErrorException find_feedback_arc_set(g, self_loops = "",
                                                      split = false)
    dfs_fas = dfs_feedback_arc_set(g)
    greedy_fas = greedy_feedback_arc_set(g)
    pagerank_fas = pagerank_feedback_arc_set(g)
    @test (1, 1) in dfs_fas && (3, 3) in dfs_fas
    @test (1, 1) in greedy_fas && (3, 3) in greedy_fas
    @test (1, 1) in pagerank_fas && (3, 3) in pagerank_fas
end
