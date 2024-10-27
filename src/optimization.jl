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

function OptProblem(graph::EdgeSubGraph)
    M = ne(graph.parent)

    l = zeros(Int, M)
    u = Int.(graph.active)
    vartypes = fill(:Bin, M)

    c = ones(Float64, M)

    return OptProblem(c, l, u, vartypes, CycleConstraint[])
end

mutable struct Solution
    status
    objval
    sol
    attrs
end

function solve_IP(O::OptProblem; solver, initial_solution = Int[],
                  use_warmstart = true, options...)
    solve_IP(Val(Symbol(solver)), O, initial_solution, use_warmstart;
             options...)
end

function solve_IP(solver::Val{T}, args...; kwargs...) where T
    error("Unknown solver `$(T)`")
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
function constrain_cycles!(O::OptProblem, cycles)
    previous_number_of_constraints = length(O.cycle_constraints)
    for cycle in cycles
        A = spzeros(Int, length(O.c))
        for edge in cycle
            A[edge] = 1
        end
        lb = 1
        ub = length(cycle)
        push!(O.cycle_constraints, CycleConstraint(A, lb, ub))
    end

    return length(O.cycle_constraints) - previous_number_of_constraints
end
