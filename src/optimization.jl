using SparseArrays: SparseArrays, SparseVector, spzeros

struct CycleConstraint
    A::SparseVector
    lb::Int
    ub::Int
end

# l and u are lower and upper bounds on the variables, i.e. the edges.
mutable struct OptProblem
    c::Vector{Float64}
    l::Vector{Int}
    u::Vector{Int}
    vartypes::Vector{Symbol}
    cycle_constraints::Vector{CycleConstraint}
end

function OptProblem(graph)
    N = nv(graph)
    M = ne(graph)

    l = zeros(Int, M)
    u = ones(Int, M)
    vartypes = fill(:Bin, M)
    edge_variables = Tuple{Int, Int}[]

    for n = 1:N
        for k in outneighbors(graph, n)
            push!(edge_variables, (n, k))
        end
    end

    c = ones(Float64, M)

    return (OptProblem(c, l, u, vartypes, CycleConstraint[]),
            edge_variables)
end

mutable struct Solution
    status
    objval
    sol
    attrs
end

function add_cycle_constraints_to_formulation(O::OptProblem)
    n = length(O.cycle_constraints)
    m = length(first(O.cycle_constraints).A)
    A = spzeros(Int, n, m)
    lb = zeros(Int, n)
    ub = zeros(Int, n)
    i = 1
    for c in O.cycle_constraints
        A[i,:] = c.A
        lb[i] = c.lb
        ub[i] = c.ub
        i += 1
    end
    return A, lb, ub
end

# Add constraints derived from cycles in the graph. The constraints
# just say that for each listed cycle, at least one edge must be
# included in the solution.
function constrain_cycles!(O::OptProblem, cycles, edges)
    previous_number_of_constraints = length(O.cycle_constraints)
    for cycle in cycles
        A = spzeros(Int, length(edges))
        cycle_edges = [(cycle[i], cycle[mod1(i + 1, length(cycle))])
                       for i in 1:length(cycle)]
        for i = 1:length(edges)
            v1, v2 = edges[i]
            if (v1, v2) in cycle_edges
                A[i] = 1
            end
        end
        lb = 1
        ub = length(cycle)
        push!(O.cycle_constraints, CycleConstraint(A, lb, ub))
    end

    return length(O.cycle_constraints) - previous_number_of_constraints
end
