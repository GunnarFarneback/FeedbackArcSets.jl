#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")

#let lemma(label_, summary, proposition, proof) = block[
  #figure(kind: "theorem", supplement: "Lemma")[] #label(label_)
  *Lemma #context counter(figure.where(kind: "theorem")).display()*, _#(summary)_: \ #proposition

  *Proof*: #proof
]

#let corollary(label_, summary, proposition, proof) = [
  #figure(kind: "theorem", supplement: "Corollary")[] #label(label_)
  *Corollary #context counter(figure.where(kind: "theorem")).display()*, _#(summary)_: \ #proposition

  *Proof*: #proof
]

#let todo(content) = [
  #box(fill: yellow, inset: 2mm)[
    #align(center)[*Under Construction*]
    #content
  ]
]

= Abstract

#pagebreak()

#outline()

#pagebreak()

= Feedback Arc Sets

= Terminology

Most of this terminology is standard but there are some
specializations for the feedback arc set problem, such as assuming
that graphs are directed or that the exact set of vertices is not
interesting in subgraphs. There are also some terms, such as partial
minimum feedback arc set, which are either invented here or uncommon
in the literature.

The exposition in this section is terse but many of the terms are
illustrated in @preliminaries.

== Quick Summary

All graphs in this document are directed and have no self-loops. A
feedback arc set is a subset of the edges of the graph, such that
removing them from the graph leaves an acyclic graph. A minimum
feedback arc set is a feedback arc set of the smallest cardinality.

== Graph

A directed graph $G = (V, E)$ is a pair of a vertex set $V$ and a
directed edge set $E subset V times V$. A directed edge $e = (v, w) in
E$ is an ordered pair of a source vertex $v in V$ and a destination
vertex $w in V$. *In this document graph is understood to mean
directed graph* and any references to undirected graphs will
explicitly use the term undirected. *The same holds for edges*. Arc is
used as a synonym for edge, but only in the context of feedback arc
sets.

== Self Loop

An edge $(v, v)$ with the same source and destination is called a self
loop. *In this document it is understood that graphs do not contain
self loops*.

== Path

A path in a graph $G = (V, E)$ is a sequence of vertices $(v_1, v_2,
..., v_n)$ such that each consecutive pair is an edge in the graph,
$(v_i, v_(i+1)) in E$ for $1 <= i < n$. A simple path is a path
without repeated vertices.

A path can equivalently be described as a sequence of connected edges.

== Cycle

A cycle in a graph $G = (V, E)$ is a path $(v_1, v_2, ..., v_n)$ where
additionally the last vertex has an edge to the first vertex, $(v_n,
v_1) in E$. A simple cycle is a cycle without repeated vertices.

A cycle can equivalently be described as a sequence of connected edges.

== Acyclic Graph

A graph $G$ is called acyclic if it does not contain any cycle. Since
graphs are understood to be directed, this is the same as a DAG, a
Directed Acyclic Graph.

== Isolated Vertex

A vertex $v in V$ is isolated in the graph $G = (V, E)$ if it is
neither the source, nor the destination, of any edge $e in E$.

== Leaf Edge

An edge $e = (v, w) in E$ is a leaf edge in a graph $G = (V, E)$ if
either the source vertex $v$ is not the destination vertex of any edge
in $E$, or conversely, the destination vertex $w$ is not the source of
any edge in $E$.

== Subgraph

In general $G' = (V', E')$ is a subgraph of $G = (V, E)$ if $G'$ is a
graph, $V' subset V$, and $E' subset E$. In this document we are
mostly concerned with edge subgraphs and unless explicitly stated
otherwise, it is assumed that $V' = V$. Subgraphs constructed by
removing edges are denoted by $G without X = (V, E without X)$, where $X subset
E$, or in the case of a single edge, $G without e = (V, E without {e})$ with $e in
E$.

== Disjoint Subgraphs

Two subgraphs $G_1 = (V, E_1)$ and $G_2 = (V, E_2)$ of $G = (V, E)$
are called disjoint if $E_1 inter E_2 = emptyset$, i.e. if they have no
edges in common.

== Feedback Arc Set <feedback-arc-set-definition>

A feedback arc set $F$ of a graph $G = (V, E)$ is a subset of the
edges, $F subset E$, such that the residual graph $G_r = G without F =
(V, E without F)$ is acyclic. I.e. removing a feedback arc set from a
graph breaks all cycles and leaves a DAG.

We use $mu(G)$ to denote the size of the smallest feedback arc set.

A feedback arc set $A$ with $|A| = mu(G)$ is called a _minimum_
feedback arc set.

A feedback arc set $A$ such that any strict subset $A' subset.neq A$
is not a feedback arc set is called a _minimal_ feedback arc set.

A partial minimum feedback arc set $P subset E$ is a subset of the
edges with the property that $mu(G) = |P| + mu(G without P)$. A subset
of a feedback arc set which forms a partial minimum feedback arc set
is also called a strong arc set.

== Strongly Connected Component

Two vertices $v_1 in V$ and $v_2 in V$ are strongly connected in a
graph $G = (V, E)$ if and only if there is both a path from $v_1$ to
$v_2$ and a path from $v_2$ to $v_1$, or equivalently, if there is a
cycle containing both vertices. This is an equivalence relation and
the equivalence classes are called strongly connected components. This
means that two vertices are strongly connected if they belong to the
same component and not strongly connected if they belong to different
components.

In the same way we can define two edges $e_1 in E$ and $e_2 in E$ to
be strongly connected if and only if there is a cycle containing both
edges. This is also an equivalence relation and we call the
equivalence classes strongly connected edge components.

= Preliminaries <preliminaries>

== Examples of Terms

#todo[Add example graphs.]

Figure 1 shows a graph $G = (V, E)$ with two vertices, $V = {v_1,
v_2}$ and two edges, $E = {e_1, e_2}$.

== Self Loops

If a graph has self loops, it depends on the application whether they
should count as a length one cycles or should be ignored because they
are irrelevant. In the former case the only way to remove the length
one cycles is to include all self loops in the feedback arc set. In
the latter case they should just be removed from the graph. Either way
they are computationally uninteresting, which is why the graphs
considered here are assumed not to have any self loops.

== Isolated Vertices

Isolated vertices cannot, obviously, be part of any cycles and thus
have no consequences for computation of feedback arc sets.

== Leaf Edges

== Feedback Arc Sets

These observations are all easy consequences of the definition of
feedback arc sets.

#lemma("lemma-subgraph-of-acyclic",
  [
  Removing edges from a graph cannot introduce new cycles
  ],
  [
    Assume that $G' = (V, E')$ is a subgraph of $G = (V, E)$ and that
    $G$ is acyclic. Then $G'$ is also acyclic.
  ],
  [
    For contradiction, assume that $G'$ is cyclic. Then there
    exists a sequence of edges $(e_1, ..., e_n)$ in $G'$ forming
    a cycle. Since $E' subset E$, this must also be a cycle in
    $G$, contradicting $G$ being acyclic.
  ]
)

#corollary("corollary-fas-in-subgraph",
  [
    A feedback arc set of a graph induces a feedback arc set of a subgraph
  ],
  [
    Assume that $G' = (V, E')$ is a subgraph of $G = (V, E)$ and that
    $F$ is a feedback arc set of $G$. Then $F inter E'$ is a feedback
    arc set of $G'$.
  ],
  [
    By the definition of a feedback arc set, $(V, E without F)$ is
    acyclic. Since $E' subset E$ it follows that
    $E' without (F inter E') subset E without F$ and by
    @lemma-subgraph-of-acyclic
    that $(V, E' without (F inter E'))$ is
    acyclic. Thus $F inter E'$ is a feedback arc set of $G'$.
  ]
)

#corollary("corollary-mu-of-subgraph",
  [
    The minimum feedback arc set cannot be larger in a subgraph
  ],
  [
    Assume that $G' = (V, E')$ is a subgraph of $G = (V, E)$. Then
    $mu(G') <= mu(G)$.
  ],
  [
    Let $F$ be a minimum feedback arc set of $G$. Then by
    @corollary-fas-in-subgraph, $F inter E'$ is a feedback arc set
    of $G'$. Since $F inter E' subset F$ it follows that
    $mu(G') <= |F inter E'| <= |F| = mu(G)$.
  ]
)

#corollary("corollary-mu-of-disjoint-subgraphs",
  [
    Splitting a graph cannot make the combined feedback arc sets larger
  ],
  [
    Assume that $G_1 = (V, E_1)$ and $G_2 = (V, E_2)$ are disjoint
    subgraphs of $G = (V, E)$, i.e. $E_1 subset E, E_2 subset E, E_1
    inter E_2 = emptyset$.
    Then $mu(G_1) + mu(G_2) <= mu(G)$.
  ],
  [

    Let $F$ be a minimum feedback arc set for $G$. By
    @corollary-fas-in-subgraph, $F_1 = E_1 inter F$ and $F_2 = E_2 inter
    F$ are feedback arc sets for $G_1$ and $G_2$ respectively. Since
    $E_1$ and $E_2$ are disjoint, $F_1$ and $F_2$ must also be
    disjoint. It follows that
    $
    mu(G_1) + mu(G_2) <= |F_1| + |F_2| <= |F| = mu(G).
    $
  ]
)

#corollary("corollary-mu-of-split-off-cycle",
  [
    Splitting off a cycle
  ],
  [
    Assume that $G' = (V, E')$ is a subgraph of $G = (V, E)$, and that
    $G'$ contains at least one cycle. Let $G_r = (V, E without E')$ be
    the residual graph.
    Then $mu(G) >= 1 + mu(G_r).$
  ],
  [
    By construction, $G'$ and $G_r$ are disjoint subgraphs of
    $G$. Since $G'$ contains at least one cycle it is not acyclic, so
    $mu(G') >= 1$. By @corollary-mu-of-disjoint-subgraphs it follows
    that
    $
      mu(G) >= mu(G') + mu(G_r) >= 1 + mu(G_r).
    $
  ]
)

#corollary("corollary-mu-of-disjoint-cycles",
  [
    Graph with $n$ disjoint cycles
  ],
  [
    Assume that $G_1, dots, G_n$ are mutually disjoint subgraphs of
    $G = (V, E)$, and that each subgraph contains at least one cycle.
    Then $mu(G) >= n$.
  ],
  [
    Repeat @corollary-mu-of-split-off-cycle $n$ times, splitting off one
    subgraph at a time.
  ]
)

#lemma("lemma-feedback-arc-in-each-cycle",
  [
    Every cycle contains a feedback arc
  ],
  [
    Assume that $F$ is a feedback arc set of $G = (V, E)$ and that $c$
    is a cycle in $G$. Then at least one of the edges in $c$ is
    contained in $F$.
  ],
  [
    For contradiction, assume that no edge in $c$ is contained in
    $F$. Then all edges in $c$ are present in the residual graph $(V,
    E without F)$, which thus contains the cycle $c$ and implies that
    $F$ is not a feedback arc set of $G$.
  ]
)

== Simple FAS bounds

These are all trivial or well-known observations.
As usual, let $G = (V, E)$ be a directed graph.

- $mu(G) <= |E|$ because $E$ is trivially a feedback arc set since the
  residual is an empty graph and thus acyclic.

- $mu(G) <= |E| - 1$. Given any $e in E$, $(E without e)$ is also
  trivially a feedback arc set since the residual contains a single
  edge and the smallest cycle requires two edges.

- $mu(G) <= (|E|) / 2$. Order $V$ arbitrarily. Let $F_f subset E$ be
  the edges for which the source vertex is ordered before the
  destination vertex (forward edges) and $F_b subset E$ be the edges
  for which the source vertex is ordered after the destination vertex
  (back edges).  These sets are complementary and it is obvious from
  the construction that both $(V, F_f)$ and $(V, F_b)$ are acyclic,
  and thus both $F_f$ and $F_b$ are feedback arc sets. Not both of
  them can be larger than $(|E|) / 2$. This bound is sharp for any
  complete directed graph, which follows from
  @corollary-mu-of-disjoint-cycles since a complete graph contains
  $(|E|) / 2$ disjoint two-cycles.

- $mu(G) >= 0$. This follows from the definition and is a sharp bound
  for acyclic graphs.

- $mu(G) >= 1$ if $G$ is not acyclic. This also follows from the
  definition since the empty set is not a feedback arc set in this
  case.

- $mu(G) >= n$ if $G$ contains $n$ two-cycles. Since two-cycles are
  necessarily disjoint, this bound follows from
  @corollary-mu-of-disjoint-cycles.

Another trivial but extremely important observation is that

- $mu(G) <= |F|$ if $F$ is a feedback arc set of $G$.

This means that upper bounds can be proved by example.

== Related Problems

There are a number of graph problems that are closely related to the
minimum feedback arc set problem or are useful to consider in this
context.

=== Minimal Feedback Arc Set

// TODO: Rewrite as an explicit algorithm?
//       Add a proof?

Finding a minimal feedback arc set (defined in
@feedback-arc-set-definition) is considerably easier than finding a
minimum feedback arc set. This can be done by starting with an
arbitrary feedback arc set, iterating over the edges, and removing any
edge that is not needed for the set to remain a feedback arc set. This
can be done in quadratic time with respect to the number of edges.

#todo[Add an example graph.]

Any minimum feedback arc set is trivially a minimal feedback arc set,
but the converse does not hold. Reducing a feedback arc set to become
minimal can be an effective way to improve an approximate solution to
the minimum feedback arc set problem.

=== Maximum Acyclic Subgraph

Finding the largest subset of edges that form a directed acyclic graph
is the equivalent complementary problem of finding a minimum feedback
arc set. I.e., if $F$ is a feedback arc set of $G = (V, E)$, then by
definition $(V, E without F)$ is acyclic and conversely if $D$ is a
subset of $E$ such that $(V, D)$ is acyclic, then $E without D$ is a
feedback arc set. Since any feedback arc set and its corresponding
maximum acyclic graph partitions $E$, it follows trivially that a
minimum $F$ corresponds to a maximum $E without F$ and vice versa.

From a practical point of view, the minimum feedback arc set is never
larger than its complement, which can be advantageous to search
algorithms.

#todo[Add an example graph for illustration.]

=== Optimal Linear Ordering

Given a linear ordering of the vertices, a feedback arc set can be
obtained as the back-edges with respect to this ordering. If this is a
minimum feedback arc set, we call the ordering optimal.

Finding an optimal ordering is equivalent to finding a minimum
feedback arc set in the sense that given an optimal ordering, a
minimum feedback arc set can immediately be derived as the back-edges,
and conversely given a minimum feedback arc set, an optimal ordering
can be found by topological ordering of the residual maximum acyclic
graph.

=== Maximum Disjoint Cycle Packing

A disjoint cycle packing is a collection of edge-disjoint cycles such
that removal of all edges in those cycles leaves an acyclic
graph. Thus it follows by definition that a disjoint cycle packing
implies a feedback arc set as the union of the edges in the cycles and
thus provides an upper bound on the size of the minimum feedback arc
set. Moreover it follows from @corollary-mu-of-disjoint-cycles that
the number of disjoint cycles provides a lower bound.

Maximizing the number of disjoint cycles is unfortunately also an
NP-hard problem and even if we can find a solution, neither the lower
nor the upper bound will in general be sharp. However, finding a
suboptimal cycle packing can be an effective way of establishing a
useful lower bound on the size of the minimum feedback arc set.

=== Weighted Minimum Feedback Arc Set

If weights are added to each edge, the minimum feedback arc set can be
generalized by searching for the feedback arc set which minimizes the
sum of the weights for the edges in the set. This is strictly more
general since setting all weights to one reduces it to the
(unweighted) minimum feedback arc set problem.

This problem is outside the scope of this document, but some of the
methods apply also for the weighted variant as is or with obvious
modifications.

=== Minimum Feedback Vertex Set

Instead of removing edges to turn a directed graph into an acyclic
graph, one can remove vertices (along with all connected edges). This
leads naturally to the minimum feedback vertex set problem, which is a
close relative of the minimum feedback arc set problem. This is also
outside of the scope of this document but some of the ideas presented
here can be expected to carry over.

= Graph Splitting

== Strongly Connected Components

== Biconnected Components

= Graph Simplification

== Two-cycles

A graph with two-cycles can always be simplified by dropping all
two-cycles. Afterwards exactly one edge from each two-cycle needs to
be added to any feedback arc set computed without the two-cycles.

#lemma("lemma-splitting-off-one-two-cycle",
  [
  Splitting off one two-cycle
  ],
  [
    Assume $G = (V, E)$ contains a two-cycle between vertices $v in V$
    and $w in V$. I.e. $e_1 = (v, w) in E$ and $e_2 = (w, v) in E$.
    Let $G_r = (V, E without {e_1, e_2})$ be the residual graph. Let
    $F_r$ be a feedback arc set for $G_r$. Then at least one of $F_r
    union e_1$ and $F_r union e_2$ is a feedback arc set for $G$.
  ],
  [
    Since $F_r$ is a feedback arc set for $G_r$, the complementary
    graph $G'_r = G_r without F_r = (V, (E without (F_r union e_1
    union e_2)))$ is acyclic. Thus $G'_r$ cannot both contain paths
    from $v$ to $w$ and from $w$ to $v$. Assume first that there is no
    path from $v$ to $w$. Then $G'_r union (w, v) = (V, (E without (F_r
    union e_1)))$ is an acyclic subgraph of $G$ and thus the
    complementary edges $F_r union e_1$ is a feedback arc set for $G$.
    Conversely, if there is no path in $G'_r$ from $w$ to $v$, then
    $F_r union e_2$ is a feedback arc set for $G$.
  ]
)

#corollary("corollary-splitting-off-two-cycles",
  [
    Splitting off multiple two-cycles
  ],
  [
    Assume $G = G_0 union G_1 union ... union G_n$, all $G_i$ are
    mutually disjoint, and $G_1, ..., G_n$ each contain a single
    two-cycle. Then \
    (1) any feedback arc set for $G_0$ can be complemented
    with one edge from each of $G_1, ..., G_n$ to form a feedback arc
    set for $G$, and \
    (2) $mu(G) = mu(G_0) + n$.
  ],
  [
    (1) follows directly by repeatedly applying
    @lemma-splitting-off-one-two-cycle. For (2) it follows from (1)
    that a feedback arc set for $G_0$ of size $mu(G_0)$ can be
    complemented by $n$ additional edges to form a feedback arc set
    for $G$ of size $mu(G_0) + n$. Thus $mu(G) <= mu(G_0) + n$.
    Conversely, repeatedly applying @corollary-mu-of-split-off-cycle
    it follows that $mu(G) >= mu(G_0) + n$.
  ]
)

This split and reconstruction can be implemented efficiently by the
following algorithm, given a graph $G = (V, E)$.

+ Remove all two-cycles from $G$ to form a reduced graph $G' = (V, E')$ and
  store the two-cycles as a collection of vertex pairs ${(v_i, w_i)}$.
+ Find a feedback arc set $F'$ for $G'$ by any method.
+ Compute a topological ordering $tau$ of the acyclic graph $(V,E'
  without F')$ so that $tau(v) < tau(w)$ for all edges $(v, w) in E'
  without F'$.
+ Initialize $F = F'$
+ For each two-cycle $(v_i, w_i)$:
  + If $tau(v_i) < tau(w_i)$, add the edge (w_i, v_i) to $F$.
  + Otherwise add the edge (v_i, w_i) to $F$.

The resulting $F$ is a feedback arc set for $G$. Moreover, if $F'$ is
a minimum feedback arc set for $G'$, then $F$ is a minimum feedback
arc set for $G$.

= Graph Decomposition

= Heuristic Methods

== Greedy Algorithm

== Pagerank Algorithm

= Exact Methods

Exact methods differ from heuristic methods in that they return proven
minimum feedback arc sets. If time limits stop them from finding an
optimal solution, they should preferably return the smallest feedback
arc set and highest lower bound found so far.

== Very Small Graphs

The feedback arc sets for graphs with at most $n$ edges can easily be
found by the method in the next section, provided that $n$ is
sufficiently small. An alternative and much faster way is to compute
feedback arc sets by static analysis of the edges, which here means
that no back-tracking search or equation solvers are involved.

In @appendix-small-graphs such methods are presented for arbitrary
graphs with a most four edges, and for single strongly connected
components with at most eight edges.

If the entire graph is very small, the value of a dedicated small
graph solver is quite limited but if graph splitting of a very large
graph produces lots of very small components, a static solver can be
helpful.

== Linear Programming Formulations

The minimum feedback arc set problem can be formulated as an integer
linear program as follows.

Let $G = (V, E)$ be the directed graph and let $CC$ be all simple
cycles in $G$. Enumerate the edges $E = {e_1, e_2, ..., e_n}$ and
associate to each edge $e_i$ a variable $x_i$, which can be
interpreted as an element of a vector $x$. This provides a natural
mapping between subsets of $E$ and binary vectors $x$, where $x_i = 1$
if and only if $e_i$ is an element in the subset.

Similarly the cycles in $CC$ can be enumerated and encoded into a
binary matrix $A$, where element $A_(i j)$ is one if and only if edge
$e_j$ is a part of cycle $i$. For later use, let $A(C)$ be the
corresponding matrix constructed from a subset of all simple cycles,
$C subset CC$.

As noted in @lemma-feedback-arc-in-each-cycle, a feedback arc set must
have at least one element in each cycle in the graph. This gives the
characterization that a binary vector $x$ corresponds to a feedback
arc set if and only if

$
sum_j A_(i j) x_j >= 1 "for all" i
$

or shorter

$
A x >= 1,
$

where the inequality is understood to hold for each element. These
constraints will be referred to as cycle constraints.

Since the size of a subset simply is the sum of the elements of the
corresponding binary vector, finding a minimum feedback arc set can
be formulated as the binary linear program
$
mu(G) =& "minimize" && sum_i x_i \
& "subject to" quad && A x >= 1, \
& &&x_i in {0, 1} "for all" i,
$ <ILP>
where any solution $x$ provides a minimum feedback arc set.

In theory this approach can find a minimum feedback arc set for any
graph but in practice it will only be effective if the graph is very
small and/or has some special structure which limits the amount of
cycles. The main problem is that the number of simple cycles tends to
increase exponentially with the size of the graph. Additionally
integer linear programming (ILP or IP) solvers tend to slow down
considerably as the problem size grows.

Even if an IP solver cannot solve to optimality, it may be able to
prove bounds on the objective function, which implies bounds on
$mu(G)$.

== Relaxations

There are two ways to relax @ILP, by considering only a subset of the
cycle constraints or by allowing non-binary $x_i$ (LP relaxization).

In the first case, let $C subset CC$ be a subset of all simple cycles
in $G$. Relax @ILP to

$
mu_(C)(G) =& "minimize" && sum_i x_i \
& "subject to" quad && A(C) x >= 1, \
& &&x_i in {0, 1} "for all" i.
$ <ILPC>

The solution $x$ still corresponds to a subset of the edges, but it is
not necessarily a feedback arc set any more. The residual graph may
contain some of the cycles in $CC without C$.

In the second case, LP relaxation, the equation becomes
$
mu^*(G) =& "minimize" && sum_i x_i \
& "subject to" quad && A x >= 1, \
& && 0 <= x_i <= 1 "for all" i.
$ <LP>

Now the solution $x$ is not necessarily possible to interpret as a
subset of the edges, since there may be fractional $x_i$.

Finally both relaxations can be applied together,

$
mu^*_(C)(G) =& "minimize" && sum_i x_i \
& "subject to" quad && A(C) x >= 1, \
& && 0 <= x_i <= 1 "for all" i.
$ <LPC>

By construction it is clear that
$
mu^*_(C)(G) &<= mu^(C)(G) <= mu(G) quad "and" \
mu^*_(C)(G) &<= mu^*(G) <= mu(G).
$

Thus all of these relaxations can be used to prove lower bounds on
$mu(G)$, regardless whether the solution $x$ provides a feedback arc
set.

For practical purposes @ILPC is the most useful relaxation, since
reducing the number of constraints allows IP solvers to work with much
larger graphs within reasonable time limits. There are two key
observations:

+ We may be lucky enough that the optimal solution $x$ provides a
  feedback arc set even if it is not guaranteed. If that is the case,
  $x$ is also a solution to @ILP.

+ If we are lucky enough, or have stumbled on a minimum feedback arc
  set by other means, $C$ only needs to contain enough cycles to make
  the lower bound sharp, in order to prove that we have a minimum
  feedback arc set, since the upper bound can always be proved by
  example. For this purpose it is also enough if the IP solver finds
  a sharp lower bound on the objective function, even if it cannot
  solve the problem to optimality.

In contrast, @LP is mostly a theoretical construction. @LPC does have
a practical value, however, because it is an LP problem with much
better complexity than the corresponding IP problem and thus can be
solved much faster. Even if the solution is fractional, @LPC provides
a lower bound and can be used for cycle generation in the next
section.

Finally it should be noted that the equality $mu^*(G) = mu(G)$ holds
for some graphs but not in general.

#todo[Add example.]

== Iterative Generation of Cycle Constraints

@ILPC can be used as the basis for a practically useful exact minimum
feedback arc set solver. The complete algorithm will be presented in
@complete-ip-algorithm. For now we focus on the structure of the
algorithm:

1. Let $C subset CC$ be an initial subset, possibly empty, of the
   simple cycles in $G$.

2. Solve @ILPC.

3. Convert the solution $x$ to an edge set $F$.

4. If $F$ is a feedback arc set, return $F$.

5. Add at least one simple cycle to $C$, which was not already in $C$.

6. Goto 2.

Subject to time and memory constraints, this is guaranteed to find an
optimal solution, since eventually $C$ will contain all of $CC$ and
when that happens, step 2 solves @ILP and step 4 returns a minimum
feedback arc set.

The hope of course, and what makes this algorithm practically useful,
is that it terminates long before $C$ encompasses all of $CC$. The key
to achieving this is the generation of cycle constraints in step 5.

=== Initial Cycles

The initial cycles can be chosen arbitrarily. Practically useful
options include:

- No cycles. Then the first iteration of step 2 can be shortcut and
  effectively the cycle generation in step 5 is used also for the
  first set of cycles.

- All cycles with a length below some threshold.

// TODO: Expand on this and refer to the section where the heuristic
// approximation is presented.
- A heuristic approximation of the disjoint cycle packing problem.

The last option is recommended but for some easily solved graphs it
might waste time.

=== Cycle Generation <cycle-generation>

Empirically new cycles provide effective constraints if they are
relatively short and not overly redundant with existing cycles.

We propose the following method to generate such cycles.

- Extend the edge set $F$ to a feedback arc set $F'$ by finding a
  feedback arc set of the residual graph $(V, E without F)$. This
  does not need to be a minimum feedback arc set and can be computed
  with any heuristic algorithm. If $F$ is non-empty, it _could_ be
  computed by a recursive call to this algorithm, but that is very
  unlikely to be an efficient use of time.

// TODO: Complete this item, check the code.
- For each edge $(v, w)$ in $F' without F$, find the shortest path
  from $w$ to $v$ which does not intersect $F$. Combined with the edge
  from $v$ to $w$, these generate

=== Time Limits

==== Non-Optimal IP Solutions

==== Lower Bounds

Every solution in step 2 provides a lower bound and the largest of
those is the final lower bound if an optimal solution is not found due
to time constraints.

==== Upper Bounds

Every feedback arc set $F'$ constructed in @cycle-generation provides
an upper bound and the smallest of those is returned as an approximate
solution if an optimal solution is not found due to time constraints.

==== Resumption

=== Complete Algorithm <complete-ip-algorithm>



= Postprocessing

= Hybrid Methods

= Benchmarks

@baharev2021
@baharev-sdopt-tearing
@wikipedia-tournament
@simpson2016
@geladaris2023
@cavallaro2023
@cavallaro2024
@dasdan-graph-benchmarks
@younger1963
@pustoslemsek

#bibliography("bibliography.yaml")

= Appendices

== Static Solvers for Small Graphs <appendix-small-graphs>

=== Up to Four Edges

=== Up to Eight Edges in a Single Strongly Connected Component

