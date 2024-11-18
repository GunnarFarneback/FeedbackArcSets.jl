import JuMP
using JuMP: MOI, @variable, @constraint, @objective

function solve_IP(::Val{:jump}, O::OptProblem, initial_solution = Int[],
                  use_warmstart = true; solver_options, solver_time, log_level)
    if !haskey(solver_options, "optimizer")
        error("The `jump` IP solver requires the backend optimizer to be passed with key `\"optimizer\"` in `solver_options`.")
    end
    optimizer = solver_options["optimizer"]
    model = JuMP.direct_model(MOI.instantiate(optimizer))
    JuMP.set_string_names_on_creation(model, false)

    # Set default solver options needed by the algorithm. These can be
    # overridden with explicit `solver_options` but that might
    # interfere with convergence or time management.
    JuMP.set_time_limit_sec(model, solver_time)
    JuMP.set_attribute(model, MOI.AbsoluteGapTolerance(), 0.0)
    JuMP.set_attribute(model, MOI.RelativeGapTolerance(), 0.0)
    log_level <= 0 && JuMP.set_silent(model)

    # Pass through solver options as raw attributes.
    for (name, value) in solver_options
        name == "optimizer" && continue
        JuMP.set_attribute(model, string(name), value)
    end

    A, lb, ub = add_cycle_constraints_to_formulation(O)
    JuMP_loadproblem!(model, A, O.l, O.u, O.c, lb, ub)

    if !isempty(initial_solution) && use_warmstart
        JuMP.set_start_value.(model[:x], initial_solution)
    end
    JuMP.optimize!(model)

    if !JuMP.is_solved_and_feasible(model)
        error("Problem not solved or feasible.")
    end

    attrs = Dict()
    attrs[:objbound] = JuMP.objective_bound(model)
    attrs[:solver] = :ip
    attrs[:model] = model
    status = JuMP.termination_status(model)
    objval = JuMP.objective_value(model)
    sol = JuMP.value.(model[:x])
    solution = Solution(status, objval, sol, attrs)
    return solution
end

function JuMP_loadproblem!(model, A, l, u, c, lb, ub)
    m, n = size(A)  # row/constraints, columns/variables
    @variable(model, x[i=1:n], Bin, lower_bound = l[i], upper_bound = u[i])
    @objective(model, Min, c' * x)
    @constraint(model, A*x .>= lb)
    @constraint(model, A*x .<= ub)

    # # If the wrapper supports Interval constraints, we can do:
    # @constraint(model, lb .<= A*x .<= ub)

    return nothing
end
