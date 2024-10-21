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
    g = copy(graph)
    while true
        old_fas_size = length(feedback_arc_set)
        for s in strongly_connected_components(g)
            length(s) == 1 && continue
            if length(s) == 2
                push!(feedback_arc_set, (s[1], s[2]))
                rem_edge!(g, s[1], s[2])
                continue
            end
            g2, I = induced_subgraph(g, s)
            (i, j) = highest_edge_pagerank(g2, num_iterations)
            e = (I[i], I[j])
            push!(feedback_arc_set, e)
            rem_edge!(g, e...)
        end
        old_fas_size == length(feedback_arc_set) && break
    end
    return feedback_arc_set
end

# The pseudo-code in the paper explicitly creates a Line Graph and
# runs page rank on that. This implementation deviates by running page
# rank directly on the edges of the original graph. Additionally it
# doesn't normalize the initial ranks to sum to one, since that
# doesn't matter for finding the maximum rank.
function highest_edge_pagerank(g, num_iterations)
    edge_ranks = Dict(e => 1.0 for e in edges(g))
    vertex_ranks = zeros(nv(g))
    for i = 1:num_iterations
        vertex_ranks .= 0
        for e in edges(g)
            vertex_ranks[e.dst] += edge_ranks[e]
        end
        for e in edges(g)
            edge_ranks[e] = vertex_ranks[e.src] / outdegree(g, e.src)
        end
    end
    e = last(findmax(last, edge_ranks))
    return e.src, e.dst
end
