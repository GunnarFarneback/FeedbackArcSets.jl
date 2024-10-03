# FeedbackArcSets.jl

`FeedbackArcSets` is a Julia package dedicated to finding the smallest
[feedback arc set](https://en.wikipedia.org/wiki/Feedback_arc_set) in
a graph. A feedback arc set is a subset of the graph's edges, such
that removing them leaves an acyclic graph.

The smallest feedback arc set problem is NP-hard, so the time needed
to find the solution grows quickly with the size of the graph, unless
it has some advantageous structure. If the problem is too big to find
an optimal solution within available time, `FeedbackArcSets` can
return the smallest found feedback arc set together with a provable
lower bound for the smallest feedback arc set.

## Exact Algorithm

    find_feedback_arc_set(graph; kwargs...)

Find the smallest feedback arc set for a directed Graphs.jl
`graph`. If time limits or other restrictions prevent finding an
optimal feedback arc set, a lower bound on the smallest feedback arc
set is returned together with the best solution found. See the
docstring for keyword arguments and return values.

## Heuristic Algorithms

These all work with a directed graph from the `Graphs` package and
return a best effort feedback arc set. The algorithms are listed in
order of increasing runtime and normally the size of the feedback arc
set decreases accordingly. All of them are much faster than finding
the exact solution.

    dfs_feedback_arc_set(graph)

Runs a DFS and returns the back edges as a feedback arc set.

    greedy_feedback_arc_set(graph; randomize = true)

Uses the greedy algorithm of Eades, Lin & Smyth [1].

    pagerank_feedback_arc_set(graph; num_iterations = 5)

Uses the Page Rank based algorithm of Geladaris, Lionakis & Tollis [2].

## Adding FeedbackArcSets

FeedbackArcSets can be added with

```
using Pkg
Pkg.add(url = "https://github.com/GunnarFarneback/FeedbackArcSets.jl.git")
```

## History

FeedbackArcSets started as a spin-off from
[LongestPaths](https://github.com/GunnarFarneback/LongestPaths.jl).

## TODO

* Make it configurable whether self-loops should be an error, ignored,
  or included in the feedback arc set.

* Register in the General registry.

* Handle non-optimal solutions from the IP solver.

* Add option to use LP instead of IP for non-optimal search.

* Add support for weighted edges.

* Generalize to other solvers than Cbc. This should be more or less
  straightforward, but is at the moment hindered by an absence of
  licenses for commercial solvers like Gurobi or CPLEX.

## References

[1] A fast and effective heuristic for the feedback arc set problem.
P Eades, X Lin, WF Smyth. Information processing letters 47 (6),
319-323, 1993.

[2] Computing a feedback arc set using pagerank.
V Geladaris, P Lionakis, IG Tollis.
International Symposium on Graph Drawing and Network Visualization,
2022. Springer
