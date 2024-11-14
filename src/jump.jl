import JuMP
using JuMP: MOI, @variable, @constraint, @objective

function solve_IP(::Val{:jump}, O::OptProblem, initial_solution = Int[],
                  use_warmstart = true; solver_options, options...)

    optimizer = pop!(solver_options, "optimizer")
    model = JuMP.direct_model(MOI.instantiate(optimizer));
    JuMP.set_string_names_on_creation(model, false)

    # Pass through solver options:
    set_JuMP_options!(model, Dict(solver_options))

    A, lb, ub = add_cycle_constraints_to_formulation(O)
    JuMP_loadproblem!(model, A, O.l, O.u, O.c, lb, ub)

    if !isempty(initial_solution) && use_warmstart
        JuMP.set_start_value.(model[:x], initial_solution)
    end
    JuMP.optimize!(model)

    if !JuMP.is_solved_and_feasible(model)
        error("Problem not solve and feasible.")
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

function set_JuMP_options!(model::JuMP.Model, options_dict::Dict)

    if haskey(options_dict, "seconds")
        JuMP.set_time_limit_sec(model, Float64(pop!(options_dict, "seconds")))
    end

    if haskey(options_dict, "allowableGap")
        JuMP.set_attribute(model, MOI.AbsoluteGapTolerance(), Float64(pop!(options_dict, "allowableGap")))
    end

    if haskey(options_dict, "relativeGap")
        JuMP.set_attribute(model, MOI.RelativeGapTolerance(), Float64(pop!(options_dict, "relativeGap")))
    end

    if haskey(options_dict, "logLevel")
        log_value = pop!(options_dict, "logLevel")
        if isinteger(log_value) && (log_value < 0 || iszero(log_value))
            @info "Solver logging set to silent."
            JuMP.set_silent(model)
        end
    end

    for (name, value) in options_dict
        JuMP.set_attribute(model, string(name), value)
    end

    # for (name, value) in options
    #     MOI.set(model, MOI.RawOptimizerAttribute(string(name)), value)
    # end

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
