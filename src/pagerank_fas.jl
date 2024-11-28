"""
    pagerank_feedback_arc_set(graph; num_iterations = 5)

Compute a feedback arc set using the Page Rank based algorithm of
Geladaris, Lionakis & Tollis. The number of page rank iterations can
be controlled by `num_iterations` and defaults to 5, which is
recommended by the paper. If the graph has self-loops, those are
included in the returned feedback arc set.

*Reference:*

Computing a feedback arc set using pagerank.
V Geladaris, P Lionakis, IG Tollis.
International Symposium on Graph Drawing and Network Visualization, 2022.
Springer
"""
function pagerank_feedback_arc_set(graph; num_iterations = 5)
    self_loops = Tuple{Int, Int}[]
    if has_self_loops(graph)
        graph = copy(graph)
        remove_self_loops!(graph, self_loops)
    end

    g = EdgeSubGraph(graph)
    feedback_arc_set = Int[]
    component_stack = Vector{Int}[]
    tarjan_workspace = TarjanWorkspace(g)
    pagerank_workspace = PagerankWorkspace(g)

    # Split the graph into strongly connected edge components. The
    # feedback arc set does not need any edges going between
    # components and each component can be analyzed in isolation.
    edge_tarjan!(component_stack, g, tarjan_workspace)

    # Analyze the components one after the other. This involves
    # finding one edge to add to the feedback arc set and remove from
    # the component. At the end the component is split into new
    # connected edge components and those are put on the stack. The
    # process ends because eventually the remainder is acyclic and
    # doesn't have any strongly connected edge components.
    while !isempty(component_stack)
        component = pop!(component_stack)
        if length(component) == 2
            # A component of size two can only be a single cycle of
            # length two and either edge must go into the feedback arc
            # set.
            push!(feedback_arc_set, first(component))
            continue
        end

        set_active_edges!(g, component)
        e = highest_edge_pagerank(g, num_iterations, pagerank_workspace)
        push!(feedback_arc_set, e)
        deactivate_edge!(g, e)
        edge_tarjan!(component_stack, g, tarjan_workspace)
    end

    if isempty(self_loops)
        return vertex_tuple.(g, feedback_arc_set)
    else
        return vcat(self_loops, vertex_tuple.(g, feedback_arc_set))
    end
end

struct PagerankWorkspace
    vertex_ranks::Vector{Float64}
    edge_ranks::Vector{Float64}
    outdegrees::Vector{Int}
end

function PagerankWorkspace(graph::EdgeSubGraph)
    return PagerankWorkspace(zeros(Float64, nv(graph.parent)),
                             zeros(Float64, ne(graph.parent)),
                             zeros(Int, nv(graph.parent)))
end

# The pseudo-code in the paper explicitly creates a Line Graph and
# runs page rank on that. This implementation deviates by running page
# rank directly on the edges of the original graph. Additionally it
# doesn't normalize the initial ranks to sum to one, since that
# doesn't matter for finding the maximum rank.
function highest_edge_pagerank(g::EdgeSubGraph, num_iterations,
                               workspace::PagerankWorkspace)
    (; vertex_ranks, edge_ranks, outdegrees) = workspace

    outdegrees .= 0
    for e in edge_list(g)
        edge_ranks[e] = 1
        outdegrees[edge_source(g, e)] += 1
    end

    for _ in 1:num_iterations
        vertex_ranks .= 0
        for e in edge_list(g)
            vertex_ranks[edge_destination(g, e)] += edge_ranks[e]
        end
        for e in edge_list(g)
            v = edge_source(g, e)
            edge_ranks[e] = vertex_ranks[v] / outdegrees[v]
        end
    end

    max_rank = 0
    best_edge = 0
    for e in edge_list(g)
        if edge_ranks[e] > max_rank
            max_rank = edge_ranks[e]
            best_edge = e
        end
    end

    return best_edge
end

struct TarjanWorkspace
    index::Vector{Int}
    lowlink::Vector{Int}
    stack::Vector{Int}
    dfs_stack::Vector{Int}
end

function TarjanWorkspace(graph::EdgeSubGraph)
    n = nv(graph.parent)
    return TarjanWorkspace(zeros(Int, n), zeros(Int, n), Int[], Int[])
end

# This is a modified version of Tarjan's algorithm to find strongly
# connected components. The basic idea is to compute strongly
# connected components of edges instead of strongly connected
# components of vertices, effectively combining strongly connected
# components and induced subgraph computations.
#
# Algorithm-wise this means three modifications of Tarjan's algorithm.
#
# 1. The outer loop is over edges instead of vertices.
#
# 2. The algorithm internally still computes strongly connected vertex
#    components but needs to output strongly connected edge
#    components. The key observation is that when a strongly connected
#    vertex component is ready to be collected, the outgoing edges
#    pointing to vertices on the stack are exactly the strongly
#    connected edge component.
#
# 3. Strongly connected components containing a single vertex
#    correspond to an empty set of edges and can be ignored.
#
# Nothing is returned from this function, components are pushed to the
# `components` vector.
function edge_tarjan!(components, g::EdgeSubGraph, workspace::TarjanWorkspace)
    (; index, lowlink, stack, dfs_stack) = workspace

    # index 0 means unvisited, index > 0 means it is on the stack, and
    # index < 0 that it has already been collected.
    # TODO: Check whether having a separate on_stack indicator is a
    # cheaper solution. Maybe this tries to be overly smart.
    index .= 0
    lowlink .= 0
    next_index = 1

    for e in edge_list(g)
        v = edge_source(g, e)
        index[v] != 0 && continue
        push!(dfs_stack, v)
        while !isempty(dfs_stack)
            v = dfs_stack[end]
            if index[v] == 0
                index[v] = next_index
                lowlink[v] = next_index
                next_index += 1
                push!(stack, v)
            end
            vertex_done = true
            for e in outedges(g, v)
                w = edge_destination(g, e)
                if index[w] == 0
                    push!(dfs_stack, w)
                    vertex_done = false
                    break
                elseif index[w] > 0
                    lowlink[v] = min(lowlink[v], abs(lowlink[w]))
                end
            end

            vertex_done || continue
            pop!(dfs_stack)
            if lowlink[v] == index[v]
                if stack[end] == v
                    # Size one component. Pop it from the stack and ignore it.
                    pop!(stack)
                    index[v] = -index[v]
                else
                    component = Int[]
                    i = 0
                    # The component collection is split into two
                    # loops, where the first collects outgoing edges
                    # pointing back to the vertex component and the
                    # second pops the stack.
                    while true
                        u = stack[end - i]
                        i += 1
                        for e in outedges(g, u)
                            w = edge_destination(g, e)
                            if index[w] > 0
                                push!(component, e)
                            end
                        end
                        u == v && break
                    end
                    while true
                        u = pop!(stack)
                        index[u] = -index[u]
                        u == v && break
                    end
                    push!(components, component)
                end
            end
        end
    end

    return
end
