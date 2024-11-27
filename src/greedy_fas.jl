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
    deltas = Dict{Int, Int}()
    prev = zeros(Int, nv(graph))
    next = zeros(Int, nv(graph))
    indegrees = indegree.(Ref(graph), vertices(graph))
    outdegrees = outdegree.(Ref(graph), vertices(graph))
    max_delta = 0
    for v in vertices(graph)
        if indegrees[v] == 0 || outdegrees[v] == 0
            delta = -2
        else
            delta = max(-1, outdegrees[v] - indegrees[v])
        end
        if delta > max_delta
            max_delta = delta
        end
        push_delta!(deltas, prev, next, delta, v)
    end

    arc_set = Tuple{Int, Int}[]

    while true
        if get(deltas, -2, 0) != 0
            class = -2
        else
            while get(deltas, max_delta, 0) == 0
                max_delta -= 1
                max_delta < 0 && break
            end
            max_delta < 0 && break
            class = max_delta
        end

        v = deltas[class]
        if randomize && class != -2
            # Note: This is far from a uniform distribution over the
            # elements in the selected class, but it doesn't have to be in
            # order to provide a meaningful randomization.
            while rand() < 0.5
                v = next[v]
                deltas[class] = v
            end
        end
        delete_delta!(deltas, prev, next, class, v)

        if outdegrees[v] > 0
            for w in outneighbors(graph, v)
                indegrees[w] == 0 && continue
                outdegrees[w] == 0 && continue
                delta = max(-1, outdegrees[w] - indegrees[w])
                indegrees[w] -= 1
                new_delta = max(-1, outdegrees[w] - indegrees[w])
                if indegrees[w] == 0
                    new_delta = -2
                end
                if new_delta != delta
                    delete_delta!(deltas, prev, next, delta, w)
                    push_delta!(deltas, prev, next, new_delta, w)
                    if new_delta > max_delta
                        max_delta = new_delta
                    end
                end
            end
            outdegrees[v] = 0
        end

        if indegrees[v] > 0
            for w in inneighbors(graph, v)
                if class != -2 && outdegrees[w] > 0
                    push!(arc_set, (w, v))
                end
                indegrees[w] == 0 && continue
                outdegrees[w] == 0 && continue
                delta = max(-1, outdegrees[w] - indegrees[w])
                outdegrees[w] -= 1
                new_delta = max(-1, outdegrees[w] - indegrees[w])
                if outdegrees[w] == 0
                    new_delta = -2
                end
                if new_delta != delta
                    delete_delta!(deltas, prev, next, delta, w)
                    push_delta!(deltas, prev, next, new_delta, w)
                    if new_delta > max_delta
                        max_delta = new_delta
                    end
                end
            end
            indegrees[v] = 0
        end
    end

    return arc_set
end

function push_delta!(deltas, prev, next, delta, v)
    if get(deltas, delta, 0) == 0
        deltas[delta] = v
        prev[v] = v
        next[v] = v
    else
        u = deltas[delta]
        w = next[u]
        prev[v] = u
        next[v] = w
        prev[w] = v
        next[u] = v
    end
end

function delete_delta!(deltas, prev, next, delta, v)
    if next[v] == v
        deltas[delta] = 0
    else
        u = prev[v]
        w = next[v]
        if deltas[delta] == v
            deltas[delta] = w
        end
        next[u] = w
        prev[w] = u
    end
end
