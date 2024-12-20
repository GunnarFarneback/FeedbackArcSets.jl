# Compute a fast feedback arc set. This doesn't have to produce a
# small feedback arc set but it must be fast and it must not include
# any edge that is not part of any cycle. As a corollary it must
# return an empty arc set for an acyclic graph.
#
# The algorithm used here is to just run DFS until all vertices are
# covered and include all found back edges into the feedback arc
# set. If the graph has self-loops, those are included in the returned
# feedback arc set.
"""
    dfs_feedback_arc_set(graph)

Compute a feedback arc set by running DFS and recording all back edges.
"""
function dfs_feedback_arc_set(graph)
    marks = zeros(nv(graph))
    feedback_arc_set = Tuple{Int, Int}[]
    for vertex in vertices(graph)
        marks[vertex] == 0 || continue
        _dfs_feedback_arc_set(graph, marks, feedback_arc_set, vertex)
    end
    return feedback_arc_set
end

# TODO: Rewrite without recursion so it doesn't stack overflow for
# deep graphs.
function _dfs_feedback_arc_set(graph, marks, feedback_arc_set, vertex)
    marks[vertex] = 1
    for neighbor in outneighbors(graph, vertex)
        if marks[neighbor] == 1
            push!(feedback_arc_set, (vertex, neighbor))
        elseif marks[neighbor] == 0
            _dfs_feedback_arc_set(graph, marks, feedback_arc_set, neighbor)
        end
    end
    marks[vertex] = 2
end

function dfs_feedback_arc_set(graph::EdgeSubGraph)
    marks = zeros(nv(graph.parent))
    feedback_arc_set = Int[]
    for edge in edge_list(graph)
        vertex = edge_source(graph, edge)
        marks[vertex] == 0 || continue
        _dfs_feedback_arc_set(graph, marks, feedback_arc_set, vertex)
    end
    return feedback_arc_set
end

function _dfs_feedback_arc_set(graph::EdgeSubGraph, marks, feedback_arc_set,
                               vertex)
    marks[vertex] = 1
    for e in outedges(graph, vertex)
        neighbor = edge_destination(graph, e)
        if marks[neighbor] == 1
            push!(feedback_arc_set, e)
        elseif marks[neighbor] == 0
            _dfs_feedback_arc_set(graph, marks, feedback_arc_set, neighbor)
        end
    end
    marks[vertex] = 2
end
