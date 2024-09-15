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

At this time two search functions are provided:

    find_feedback_arc_set(graph; kwargs...)

Find the smallest feedback arc set for a directed Graphs.jl
`graph`. If time limits or other restrictions prevent finding an
optimal feedback arc set, a lower bound on the smallest feedback arc
set is returned together with the best solution found. See the
docstring for keyword arguments and return values.

    fast_feedback_arc_set(graph)

Returns a non-optimal (and probably rather far from optimum) feedback
arc set for a directed Graphs.jl `graph`. This is primarily intended
for internal use but could be useful if you don't need a very good
solution.

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

* Add tests.

* Make it configurable whether self-loops should be an error, ignored,
  or included in the feedback arc set.

* Register in the General registry.

* Handle non-optimal solutions from the IP solver.

* Add option to use LP instead of IP for non-optimal search.

* Pick up cycles to add as constraints in a smarter way.

* Add support for weighted edges.

* Generalize to other solvers than Cbc. This should be more or less
  straightforward, but is at the moment hindered by an absence of
  licenses for commercial solvers like Gurobi or CPLEX.
