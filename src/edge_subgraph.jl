# This is a graph representation in terms of static parent graph and a
# dynamic subset of edges from the parent graph. The vertices of the
# subgraph are induced from what the edges touch. Vertices and edges
# are always referred to by their numbers in the parent graph.
mutable struct EdgeSubGraph
    # Parent graph. Numbering of vertices and edges in the subgraph
    # always relate to the parent graph.
    parent::SimpleDiGraph

    # True if and only if the corresponding edge is included in the
    # subgraph.
    active::Vector{Bool}

    # List of edges in the subgraph. Elements must be unique but may
    # come in any order. This information is redundant with `active`
    # and the reason both are present is that `edge_list` is efficient
    # for iteration whereas `active` is efficient for membership
    # testing.
    edge_list::Vector{Int}

    # Source vertex of each edge in the parent graph.
    src::Vector{Int}

    # Destination vertex of each edge in the parent graph.
    dst::Vector{Int}

    # List of outgoing edges for each vertex in the parent graph. To
    # loop over outgoing edges within the subgraph, the `active` field
    # must be consulted. The `outedges` iterator does this for you.
    fadj::Vector{Vector{Int}}
end

# Let the graph act as a scalar in broadcasting.
Base.Broadcast.broadcastable(graph::EdgeSubGraph) = Ref(graph)

# Constructor from a SimpleDiGraph and an optional edge_list. By
# default all edges are included in the subgraph.
function EdgeSubGraph(graph::SimpleDiGraph, edge_list = 1:ne(graph))
    active = fill(false, ne(graph))
    for e in edge_list
        active[e] = true
    end
    src = [e.src for e in edges(graph)]
    dst = [e.dst for e in edges(graph)]
    fadj = [Int[] for _ in 1:nv(graph)]
    for (i, e) in enumerate(edges(graph))
        push!(fadj[e.src], i)
    end
    return EdgeSubGraph(graph, active, edge_list, src, dst, fadj)
end

# Replace the edges of the subgraph with a new subset of the parent
# graph's edges.
function set_active_edges!(graph::EdgeSubGraph, edge_list::Vector{Int})
    graph.edge_list = edge_list
    graph.active .= false
    for e in edge_list
        graph.active[e] = true
    end
    return
end

# Remove one edge from the subgraph.
function deactivate_edge!(graph::EdgeSubGraph, edge::Int)
    i = findfirst(==(edge), graph.edge_list)
    isnothing(i) && return
    if i < length(graph.edge_list)
        graph.edge_list[i] = pop!(graph.edge_list)
    else
        pop!(graph.edge_list)
    end
    graph.active[edge] = false
    return
end

# Remove multiple edges from the subgraph.
function deactivate_edges!(graph::EdgeSubGraph, edges::Vector{Int})
    for edge in edges
        deactive_edge!(graph, edge)
    end
    return
end

# Accessors for struct fields.
edge_list(graph::EdgeSubGraph) = graph.edge_list
edge_is_active(graph::EdgeSubGraph, edge::Int) = graph.active[edge]
edge_source(graph::EdgeSubGraph, edge::Int) = graph.src[edge]
edge_destination(graph::EdgeSubGraph, edge::Int) = graph.dst[edge]
vertex_tuple(graph::EdgeSubGraph, edge::Int) = (graph.src[edge],
                                                graph.dst[edge])

# Iterator of subgraph edges going out from a given vertex.
#
# Note: It's slightly faster do this yourself as
#
#     for e in graph.fadj[vertex]
#         graph.active[vertex] || continue
#
# instead of
#
#     for e in outedges(graph, vertex)
#
# but is obviously less nice.
outedges(graph::EdgeSubGraph, vertex::Int) = OutEdges(graph.fadj[vertex],
                                                      graph.active)

struct OutEdges
    fadj::Vector{Int}
    active::Vector{Bool}
end

@inline function Base.iterate((; fadj, active)::OutEdges, i = 1)
    while i <= length(fadj)
        @inbounds e = fadj[i]
        i += 1
        if active[e]
            return (e, i)
        end
    end
    return nothing
end

Base.IteratorSize(::Type{OutEdges}) = Base.SizeUnknown()
