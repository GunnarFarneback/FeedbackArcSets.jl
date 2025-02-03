"""
    FeedbackArcSets

Find the smallest feedback arc set in a directed graph. The smallest
feedback arc set problem is NP-hard, so the time needed to find the
solution grows quickly with the size of the graph, unless it has some
advantageous structure.

    find_feedback_arc_set(graph; kwargs)

Find the smallest feedback arc set in a Graphs directed `graph`.
"""
module FeedbackArcSets

export FeedbackArcSet, find_feedback_arc_set, dfs_feedback_arc_set,
       greedy_feedback_arc_set, pagerank_feedback_arc_set,
       is_feedback_arc_set,
       baharev_benchmark, dasdan_benchmark, snap_benchmark,
       small_graphs_benchmark

using Graphs: Graphs, SimpleDiGraph, add_edge!, edges, has_edge, has_self_loops,
              ne, nv, outneighbors, rem_edge!, simplecycles_iter,
              simplecycles_limited_length, vertices, a_star,
              inneighbors, outneighbors, indegree, outdegree,
              strongly_connected_components, induced_subgraph, is_cyclic
using Printf: Printf, @printf

include("edge_subgraph.jl")
include("optimization.jl")
include("cbc.jl")
include("highs.jl")
include("jump.jl")
include("dfs_fas.jl")
include("greedy_fas.jl")
include("pagerank_fas.jl")
include("benchmarks.jl")

"""
    FeedbackArcSet

Type used for the return values of `feedback_arc_set`. See the function
documentation for more information.
"""
mutable struct FeedbackArcSet
    lower_bound::Int
    feedback_arc_set::Vector{Tuple{Int, Int}}
    cycles::Vector{Vector{Int}}
end

function Base.show(io::IO, x::FeedbackArcSet)
    upper_bound = length(x.feedback_arc_set)
    if x.lower_bound == upper_bound
        println(io, "Optimal Feedback arc set of size $(x.lower_bound).")
    else
        println(io, "Feedback arc set with lower bound $(x.lower_bound) and upper bound $(upper_bound).")
    end
end

# Corresponding struct used internally, returned from
# `find_feedback_arc_set_ip`. The difference is that edges are
# represented by integers instead of two-tuples and that cycles are
# represented by edges instead of vertices.
mutable struct InternalFeedbackArcSet
    lower_bound::Int
    feedback_arc_set::Vector{Int}
    cycles::Vector{Vector{Int}}
end

"""
    find_feedback_arc_set(graph)

Find the smallest feedback arc set for `graph`, which must be a
directed graph from the `Graphs` package.

    find_feedback_arc_set(graph; kwargs...)

By adding keyword arguments it is possible to guide the search or
obtain non-optimal solutions and bounds in shorter time than the full
solution.

*Keyword arguments:*

* `self_loops`: What to do with self-loops in the graph. The options
  are `"error"` to give an error, `"ignore"` to ignore them, or
  `"include"` to include them in the feedback arc set. The default is
  to give an error.

* `split`: Whether to apply graph splitting. If true the graph is
  split into strongly connected components, which are analyzed
  separately. Default is true but you might want to set it to false if
  you know that the graph only consists of one component or to
  simplify time management. Custom iteration callbacks may also be
  simpler to manage without splitting.

* `max_iterations`: Stop search after this number of iterations. If
  graph splitting is activated, the number of iterations applies to
  each component. Defaults to a very high number.

* `time_limit`: Stop search after this number of seconds. If graph
  splitting is activated the remaining time is split proportionally to
  the number of edges in the components, with a minimum of one
  second. Defaults to a very high number.

* `solver`: IP solver. Currently `"cbc"` (default) and `"highs"` are
  supported.

* `solver_options`: Dictionary of options passed to the IP solver.
  These are specific to the solver and override settings used by the
  algorithm. In particular setting gaps can interfere with the
  convergence of the search in non-obvious ways and setting time
  limits can interfere with the time management.
  When using the `jump` solver, the backend optimizer must be passed
  with the `"optimizer"` key.

* `solver_time_limit`: Maximum time to spend in each iteration to
  solve integer programs, measured in seconds. This will be gradually
  increased if the IP solver does not find a useful solution in the
  allowed time. Defaults to 10 seconds.

* `initial_arc_set`: Starting arc set as a vector of edges, where
  each edge is represented by a tuple of two `Int`s. This does not
  have to be a complete feedback arc set. This is used to guide the
  solution and might improve the solution speed. Edges which are not
  present in the graph are silently ignored. Default is an empty
  vector.

* `initial_cycles`: Initial set of cycles to use as contraints by the
  IP solver. Each cycle is given as a vector of vertices. Invalid
  cycles are silently dropped. Default is no initial cycles.

* `log_level`: Amount of verbosity during search. 0 is maximally
  quiet. Level 1 prints results of static analysis and if there are
  multiple searches a summary of those. Level 2 adds printing of
  progress for each search iteration. Higher values add diagnostics
  from the IP solver calls. Default value is 2.

* `iteration_callback`: A function provided here is called during each
  iteration. The function should take one argument which is a named
  tuple with diagnostic information, see the code for exact
  specifications. Return `true` to continue search and `false` to stop
  search. The default is the `print_iteration_data` function.

The return value is of the `FeedbackArcSet` type and contains the
following fields:

* `lower_bound`: Lower bound for feedback arc set.

* `feedback_arc_set`: Vector of the edges in the smallest found
  feedback arc set. Each edge is represented as a tuple of two `Int`s.

* `internals`: Dict containing a variety of information about the
  search. If graph splitting is active, this information only applies
  to the last searched component.

"""
function find_feedback_arc_set(graph;
                               self_loops = "error",
                               split = true,
                               time_limit = typemax(Int),
                               initial_arc_set = Tuple{Int, Int}[],
                               initial_cycles = Vector{Int}[],
                               log_level = 2,
                               kwargs...)
    self_loops ∈ ["error", "ignore", "include"] ||
        error("`self_loops` must be one of \"error\", \"ignore\", \"include\".")

    # Remove self loops if there are any.
    removed_self_loops = Tuple{Int, Int}[]
    if has_self_loops(graph)
        self_loops == "error" &&
            error("Self-loops found in the graph. Remove them from the graph or use the `self_loops` keyword argument.")
        graph = copy(graph)
        remove_self_loops!(graph, removed_self_loops)
    end

    edge_graph = EdgeSubGraph(graph)
    solution = FeedbackArcSet(0, Tuple{Int, Int}[], Vector{Int}[])

    if split == false
        initial_edges = initial_arc_set_to_edges(edge_graph, initial_arc_set)
        initial_edge_cycles =
            vertex_cycles_to_edge_cycles(edge_graph, initial_cycles)
        solution_ = find_feedback_arc_set_ip(edge_graph; time_limit,
                                             initial_solution = initial_edges,
                                             initial_cycles = initial_edge_cycles,
                                             log_level, kwargs...)
        solution.lower_bound = solution_.lower_bound
        solution.feedback_arc_set = vertex_tuple.(Ref(edge_graph),
                                                  solution_.feedback_arc_set)
        solution.cycles = edge_cycles_to_vertex_cycles(edge_graph,
                                                       solution_.cycles)
    else
        start_time = time()
        components = Vector{Int}[]
        tarjan_workspace = TarjanWorkspace(edge_graph)
        edge_tarjan!(components, edge_graph, tarjan_workspace)
        sort!(components, by = length)
        component_sizes = length.(components)
        if log_level >= 1
            for i in unique(component_sizes)
                n = count(==(i), component_sizes)
                println("Found $n component(s) of size $i")
            end
        end
        lower_bound = 0
        feedback_arcs = Tuple{Int, Int}[]
        for (i, component) in enumerate(components)
            set_active_edges!(edge_graph, component)
            remaining_time = time_limit - (time()- start_time)
            remaining_edges = sum(@view component_sizes[i:end])
            work_fraction = component_sizes[i] / remaining_edges
            time_quota = max(1, remaining_time * work_fraction)
            if log_level >= 1
                println("  Analyzing component $i of size $(component_sizes[i]):")
            end
            initial_edges = initial_arc_set_to_edges(edge_graph,
                                                     initial_arc_set)
            initial_edge_cycles =
                vertex_cycles_to_edge_cycles(edge_graph, initial_cycles)
            component_solution =
                find_feedback_arc_set_ip(edge_graph;
                                         time_limit = time_quota,
                                         initial_solution = initial_edges,
                                         initial_cycles = initial_edge_cycles,
                                         log_level, kwargs...)
            if log_level >= 1
                lb = component_solution.lower_bound
                ub = length(component_solution.feedback_arc_set)
                println("    Bounds: [$(lb), $(ub)]")
            end
            lower_bound += component_solution.lower_bound
            append!(feedback_arcs, vertex_tuple.(Ref(edge_graph),
                                                 component_solution.feedback_arc_set))
            append!(solution.cycles,
                    edge_cycles_to_vertex_cycles(edge_graph,
                                                 component_solution.cycles))
        end
        solution.lower_bound = lower_bound
        solution.feedback_arc_set = feedback_arcs
    end

    if self_loops == "include"
        append!(solution.feedback_arc_set, removed_self_loops)
        solution.lower_bound += length(removed_self_loops)
    end

    return solution
end

# Convert edges represented as tuples to the corresponding Int in the
# edge graph. Additionally this drops any edges which are not present
# and active in the edge graph.
function initial_arc_set_to_edges(edge_graph, initial_arc_set)
    initial_edges = Int[]
    for (v, w) in initial_arc_set
        for e in outedges(edge_graph, v)
            if edge_destination(edge_graph, e) == w
                push!(initial_edges, e)
                break
            end
        end
    end
    return initial_edges
end

function find_feedback_arc_set_ip(graph::EdgeSubGraph;
                                  max_iterations = typemax(Int),
                                  time_limit = typemax(Int),
                                  solver = "cbc",
                                  solver_options = Dict{String, Any}(),
                                  solver_time_limit = 10,
                                  log_level = 2,
                                  initial_solution = Int[],
                                  initial_cycles = Vector{Int}[],
                                  iteration_callback = print_iteration_data)
    O = OptProblem(graph)

    cycles = short_cycles_through_given_edges(graph,
                                              dfs_feedback_arc_set(graph))

    # If no cycles found, return the trivial optimum.
    if isempty(cycles)
        return InternalFeedbackArcSet(0, Tuple{Int, Int}[], Vector{Int}[])
    end

    if isempty(initial_cycles)
        constrain_cycles!(O, cycles)
    else
        constrain_cycles!(O, initial_cycles)
    end

    lower_bound = 1
    solution = nothing
    
    start_time = time()
    best_arc_set = complement_initial_solution!(graph, initial_solution)
    local arc_set
    
    for iteration = 1:max_iterations
        solver_time = min(solver_time_limit, time_limit - (time() - start_time))
        if (solver_time < solver_time_limit / 2
            && iteration > 1
            && time() - start_time > 0.5 * time_limit)

            break
        end

        best_solution = arcs_to_solution(graph, best_arc_set)
        solution = solve_IP(O; solver, solver_options,
                            initial_solution = best_solution, solver_time,
                            log_level = log_level - 2)
        
        objbound = solution.attrs[:objbound]

        arc_set, cycles = extract_arc_set_and_cycles(graph, solution.sol)

        if length(arc_set) < length(best_arc_set)
            best_arc_set = arc_set
        end

        lower_bound = max(lower_bound, objbound)
        rounded_lower_bound = ceil(Int, lower_bound - 0.01)

        iteration_data = (log_level = log_level,
                          iteration = iteration,
                          elapsed_time = time() - start_time,
                          lower_bound = lower_bound,
                          solver_time = solver_time,
                          num_cycles = length(O.cycle_constraints),
                          solution = solution,
                          objbound = objbound,
                          arc_set = arc_set,
                          best_arc_set = best_arc_set,
                          cycles = cycles)

        if !iteration_callback(iteration_data)
            break
        end

        if iteration == max_iterations
            break
        end

        if length(best_arc_set) == rounded_lower_bound
            # Optimum found.
            break
        end

        if length(cycles) == 0
            # This can only happen if IP solver did not solve to
            # optimality. Increase the solver time.
            solver_time_limit *= 1.6
        end

        constrain_cycles!(O, cycles)
    end

    return InternalFeedbackArcSet(ceil(Int, lower_bound - 0.01),
                                  best_arc_set, extract_cycles(O))
end

function print_iteration_data(data)
    if data.log_level >= 2
        @printf("%3d %5d ", data.iteration, round(Int, data.elapsed_time))
        print("[$(data.lower_bound) $(length(data.best_arc_set))] ")
        println("$(data.num_cycles) $(data.solution.status) $(data.solution.objval) $(data.objbound)")
    end
    return true
end

function arcs_to_solution(graph, arcs)
    solution = zeros(Int, ne(graph.parent))
    for arc in arcs
        solution[arc] = 1
    end
    return solution
end

# The initial solution may not be a complete feedback arc set,
# especially the default empty set. This function makes it complete.
function complement_initial_solution!(graph, arc_set)
    saved_edges = copy(edge_list(graph))
    for edge in arc_set
        deactivate_edge!(graph, edge)
    end
    append!(arc_set, dfs_feedback_arc_set(graph))
    set_active_edges!(graph, saved_edges)
    return arc_set
end


# Extract the possibly partial feedback arc set from the solution and
# find some additional cycles if it is incomplete.
function extract_arc_set_and_cycles(graph, solution)
    saved_edges = copy(edge_list(graph))
    arc_set = []
    for e in 1:length(solution)
        if solution[e] >= 0.5
            push!(arc_set, e)
            deactivate_edge!(graph, e)
        end
    end

    additional_arcs = dfs_feedback_arc_set(graph)
    append!(arc_set, additional_arcs)
    cycles = short_cycles_through_given_edges(graph, additional_arcs)

    set_active_edges!(graph, saved_edges)

    return arc_set, cycles
end

# Given that (v, w) is part of at least one cycle, we can find the
# shortest cycle it is part of by computing the shortest path from w
# to v.
function short_cycles_through_given_edges(graph, edges)
    return [short_cycle_through_given_edge(graph, e) for e in edges]
end

function short_cycle_through_given_edge(graph, edge)
    v = edge_source(graph, edge)
    w = edge_destination(graph, edge)
    cycle = [edge]
    queue = [w]
    edges = zeros(Int, nv(graph.parent))
    offset = edge
    while !isempty(queue)
        u = popfirst!(queue)
        n = length(graph.fadj[u])
        for i in 1:n
            # There are usually multiple paths of the same length. If
            # we always add neighbors in the same order, some edges
            # are more likely to be part of the returned cycle than
            # others. This is problematic since more diversity
            # improves the IP lower bound when the cycles are used as
            # constraints. We obtain this diversity by offsetting the
            # search order differently for each starting edge.
            e = graph.fadj[u][mod1(i + offset, n)]
            edge_is_active(graph, e) || continue
            t = edge_destination(graph, e)
            if t == v
                edges[t] = e
                empty!(queue)
                break
            elseif edges[t] == 0
                edges[t] = e
                push!(queue, t)
            end
        end
    end
    e = edges[v]
    while true
        pushfirst!(cycle, e)
        edge_source(graph, e) == w && break
        e = edges[edge_source(graph, e)]
    end
    return cycle
end

"""
    is_feedback_arc_set(graph, arc_set)

Test whether a set of edges, given as an iterable of two-element
tuples, is a feedback arc set. This means that removing the edges in
`arc_set` from `graph` leaves an acyclic graph.
"""
function is_feedback_arc_set(graph, arc_set)
    graph = copy(graph)
    for (v, w) in arc_set
        rem_edge!(graph, v, w)
    end
    return !is_cyclic(graph)
end

function remove_self_loops!(graph, self_loops)
    for v in vertices(graph)
        if has_edge(graph, v, v)
            rem_edge!(graph, v, v)
            push!(self_loops, (v, v))
        end
    end
    return
end

end
