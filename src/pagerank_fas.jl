# TODO: Investigate performance and optimize.
"""
    pagerank_feedback_arc_set(graph; num_iterations = 5)

Compute a feedback arc set using the Page Rank based algorithm of
Geladaris, Lionakis & Tollis. The number of page rank iterations can
be controlled by `num_iterations` and defaults to 5, which is
recommended by the paper.

*Reference:*

Computing a feedback arc set using pagerank.
V Geladaris, P Lionakis, IG Tollis.
International Symposium on Graph Drawing and Network Visualization, 2022.
Springer

"""
function pagerank_feedback_arc_set(graph; num_iterations = 5)
    feedback_arc_set = Tuple{Int, Int}[]

    # TODO: Package up the various arrays into structs.
    edge_src = [e.src for e in edges(graph)]
    edge_dst = [e.dst for e in edges(graph)]
    edge_fadj = [Int[] for _ in 1:nv(graph)]
    for (i, e) in enumerate(edges(graph))
        push!(edge_fadj[e.src], i)
    end
    edges_list = collect(1:ne(graph))
    vertex_ranks = zeros(nv(graph))
    edge_ranks = zeros(ne(graph))
    outdegrees = zeros(Int, nv(graph))

    component_stack = Vector{Int}[]
    edges_ind = fill(false, ne(graph))

    index = zeros(Int, nv(graph))
    lowlink = zeros(Int, nv(graph))
    tarjan_stack = Int[]
    dfs_stack = Int[]

    edge_tarjan!(component_stack, edge_src, edge_dst, edge_fadj, edges_list,
                 edges_ind, index, lowlink, tarjan_stack, dfs_stack)

    while !isempty(component_stack)
        s = pop!(component_stack)
        if length(s) == 2
            push!(feedback_arc_set, (edge_src[s[1]], edge_dst[s[1]]))
            continue
        end

        i = highest_edge_pagerank(edge_src, edge_dst, s,
                                  vertex_ranks, edge_ranks, outdegrees,
                                  num_iterations)
        e = (edge_src[i], edge_dst[i])
        push!(feedback_arc_set, e)
        edge_tarjan!(component_stack, edge_src, edge_dst, edge_fadj,
                     s, edges_ind, index, lowlink, tarjan_stack, dfs_stack, i)
    end
    return feedback_arc_set
end

# The pseudo-code in the paper explicitly creates a Line Graph and
# runs page rank on that. This implementation deviates by running page
# rank directly on the edges of the original graph. Additionally it
# doesn't normalize the initial ranks to sum to one, since that
# doesn't matter for finding the maximum rank.
function highest_edge_pagerank(edge_src, edge_dst, edges_list,
                               vertex_ranks, edge_ranks, outdegrees,
                               num_iterations)
    outdegrees .= 0
    for e in edges_list
        edge_ranks[e] = 1
        outdegrees[edge_src[e]] += 1
    end

    for _ in 1:num_iterations
        vertex_ranks .= 0
        for e in edges_list
            vertex_ranks[edge_dst[e]] += edge_ranks[e]
        end
        for e in edges_list
            edge_ranks[e] = vertex_ranks[edge_src[e]] / outdegrees[edge_src[e]]
        end
    end
    max_rank = 0
    best_edge = 0
    for e in edges_list
        if edge_ranks[e] > max_rank
            max_rank = edge_ranks[e]
            best_edge = e
        end
    end
    return best_edge
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
function edge_tarjan!(components, edge_src, edge_dst, edge_fadj,
                      edges_list, edges_ind, index, lowlink, stack,
                      dfs_stack, forbidden = 0)
    edges_ind .= false
    for e in edges_list
        edges_ind[e] = true
    end
    if forbidden != 0
        edges_ind[forbidden] = false
    end

    # index 0 means unvisited, index > 0 means it is on the stack, and
    # index < 0 that it has already been collected.
    # TODO: Check whether having a separate on_stack indicator is a
    # cheaper solution. Maybe this tries to be overly smart.
    index .= 0
    lowlink .= 0
    next_index = 1

    for e in edges_list
        v = edge_src[e]
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
            for e in edge_fadj[v]
                edges_ind[e] || continue
                w = edge_dst[e]
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
                        for e in edge_fadj[u]
                            if edges_ind[e]
                                w = edge_dst[e]
                                if index[w] > 0
                                    push!(component, e)
                                end
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
