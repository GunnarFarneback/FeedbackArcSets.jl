# Note: This implementation can be improved with more efficient data
# structures.
"""
    greedy_feedback_arc_set(graph; randomize = true)

Compute a feedback arc set using the greedy algorithm of Eades, Lin &
Smyth. This implementation defaults to randomizing the choice of
equivalently attractive edges. Setting `randomize = false` makes it
deterministic.

*Reference:*

A fast and effective heuristic for the feedback arc set problem.
P Eades, X Lin, WF Smyth. Information processing letters 47 (6),
319-323, 1993.
"""
function greedy_feedback_arc_set(graph; randomize = true)
    head = Int[]
    tail = Int[]
    sources = Int[]
    sinks = Int[]
    deltas = Dict{Int, Set{Int}}()
    indegrees = indegree.(Ref(graph), vertices(graph))
    outdegrees = outdegree.(Ref(graph), vertices(graph))
    for v in vertices(graph)
        if indegrees[v] == 0
            push!(sources, v)
        elseif outdegrees[v] == 0
            push!(sinks, v)
        else
            push!(get!(Set{Int}, deltas, outdegrees[v] - indegrees[v]), v)
        end
    end

    nonempty_deltas = Set{Int}(keys(deltas))

    while true
        while !isempty(sources)
            v = pop!(sources)
            push!(head, v)
            for w in outneighbors(graph, v)
                indegrees[w] == 0 && continue
                outdegrees[w] == 0 && continue
                delta = outdegrees[w] - indegrees[w]
                delete!(deltas[delta], w)
                isempty(deltas[delta]) && delete!(nonempty_deltas, delta)
                indegrees[w] -= 1
                if indegrees[w] == 0
                    push!(sources, w)
                else
                    push!(get!(Set{Int}, deltas, delta + 1), w)
                    push!(nonempty_deltas, delta + 1)
                end
            end
            outdegrees[v] = 0
        end
        while !isempty(sinks)
            v = pop!(sinks)
            indegrees[v] == 0 && continue
            pushfirst!(tail, v)
            for w in inneighbors(graph, v)
                indegrees[w] == 0 && continue
                outdegrees[w] == 0 && continue
                delta = outdegrees[w] - indegrees[w]
                delete!(deltas[delta], w)
                isempty(deltas[delta]) && delete!(nonempty_deltas, delta)
                outdegrees[w] -= 1
                if outdegrees[w] == 0
                    push!(sinks, w)
                else
                    push!(get!(Set{Int}, deltas, delta - 1), w)
                    push!(nonempty_deltas, delta - 1)
                end
            end
            indegrees[v] = 0
        end
        isempty(nonempty_deltas) && break
        max_delta = maximum(nonempty_deltas)
        if randomize
            v = rand(deltas[max_delta])
            delete!(deltas[max_delta], v)
        else
            v = pop!(deltas[max_delta])
        end
        push!(head, v)
        isempty(deltas[max_delta]) && delete!(nonempty_deltas, max_delta)
        for w in outneighbors(graph, v)
            indegrees[w] == 0 && continue
            outdegrees[w] == 0 && continue
            delta = outdegrees[w] - indegrees[w]
            delete!(deltas[delta], w)
            isempty(deltas[delta]) && delete!(nonempty_deltas, delta)
            indegrees[w] -= 1
            if indegrees[w] == 0
                push!(sources, w)
            else
                push!(get!(Set{Int}, deltas, delta + 1), w)
                push!(nonempty_deltas, delta + 1)
            end
        end
        outdegrees[v] = 0
        for w in inneighbors(graph, v)
            indegrees[w] == 0 && continue
            outdegrees[w] == 0 && continue
            delta = outdegrees[w] - indegrees[w]
            delete!(deltas[delta], w)
            isempty(deltas[delta]) && delete!(nonempty_deltas, delta)
            outdegrees[w] -= 1
            if outdegrees[w] == 0
                push!(sinks, w)
            else
                push!(get!(Set{Int}, deltas, delta - 1), w)
                push!(nonempty_deltas, delta - 1)
            end
        end
        indegrees[v] = 0
    end
    pos = Dict((v => i) for (i, v) in enumerate(vcat(head, tail)))
    return [(e.src, e.dst) for e in edges(graph) if (pos[e.src] > pos[e.dst])]
end
