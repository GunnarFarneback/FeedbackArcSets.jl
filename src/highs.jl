import HiGHS
using HiGHS: Highs_create, Highs_setStringOptionValue, Highs_setSolution,
             Highs_run, Highs_getDoubleInfoValue, Highs_getSolution,
             Highs_getObjectiveValue, Highs_destroy, Highs_passMip,
             Highs_getModelStatus, HighsInt

function solve_IP(::Val{:highs}, O::OptProblem, initial_solution = Int[],
                  use_warmstart = true; solver_options, solver_time, log_level)
    model = Highs_create()
    Highs_setStringOptionValue(model, "time_limit", string(solver_time))
    Highs_setStringOptionValue(model, "mip_abs_gap", "0")
    Highs_setStringOptionValue(model, "mip_rel_gap", "0")
    Highs_setStringOptionValue(model, "log_to_console", string(log_level > 0))

    for (name, value) in solver_options
        Highs_setStringOptionValue(model, string(name), string(value))
    end

    A, lb, ub = add_cycle_constraints_to_formulation(O)
    highs_loadproblem!(model, A, O.l, O.u, O.c, lb, ub)

    if !isempty(initial_solution) && use_warmstart
        Highs_setSolution(model, initial_solution, C_NULL, C_NULL, C_NULL)
    end
    Highs_run(model)

    attrs = Dict()
    p = Ref{Cdouble}()
    Highs_getDoubleInfoValue(model, "mip_dual_bound", p)
    attrs[:objbound] = p[]
    attrs[:solver] = :ip
    sol = zeros(size(A, 2))
    Highs_getSolution(model, sol, C_NULL, C_NULL, C_NULL)
    solution = Solution(highs_status(model), Highs_getObjectiveValue(model),
                        sol, attrs)
    Highs_destroy(model)
    return solution
end

function highs_loadproblem!(model, A, l, u, c, lb, ub)
    Highs_passMip(model,                            # highs
                  size(A, 2),                       # num_col
                  size(A, 1),                       # num_row
                  length(A.nzval),                  # num_nz
                  HiGHS.kHighsMatrixFormatColwise,  # a_format
                  HiGHS.kHighsObjSenseMinimize,     # sense
                  0.0,                              # offset
                  Float64.(c),                      # col_cost
                  Float64.(l),                      # col_lower
                  Float64.(u),                      # col_upper
                  Float64.(lb),                     # row_lower
                  Float64.(ub),                     # row_upper
                  HighsInt.(A.colptr .- 1),         # a_start
                  HighsInt.(A.rowval .- 1),         # a_index
                  Float64.(A.nzval),                # a_value
                  fill(HiGHS.kHighsVarTypeInteger, size(A, 2)))  # integrality
end

function highs_status(model)
    status = Highs_getModelStatus(model)
    status == HiGHS.kHighsModelStatusOptimal && return :Optimal
    status == HiGHS.kHighsModelStatusTimeLimit && return :UserLimit
    status == HiGHS.kHighsModelStatusInfeasible && return :Infeasible
    status == HiGHS.kHighsModelStatusInterrupt && return :Error
    status == HiGHS.kHighsModelStatusIterationLimit && return :UserLimit
    status == HiGHS.kHighsModelStatusLoadError && return :Error
    status == HiGHS.kHighsModelStatusModelEmpty && return :Error
    status == HiGHS.kHighsModelStatusModelError && return :Error
    status == HiGHS.kHighsModelStatusNotset && return :InternalError
    status == HiGHS.kHighsModelStatusObjectiveBound && return :UserLimit
    status == HiGHS.kHighsModelStatusObjectiveTarget && return :UserLimit
    status == HiGHS.kHighsModelStatusPostsolveError && return :Error
    status == HiGHS.kHighsModelStatusPresolveError && return :Error
    status == HiGHS.kHighsModelStatusSolutionLimit && return :UserLimit
    status == HiGHS.kHighsModelStatusSolveError && return :Error
    status == HiGHS.kHighsModelStatusUnbounded && return :InternalError
    status == HiGHS.kHighsModelStatusUnboundedOrInfeasible && return :InternalError
    status == HiGHS.kHighsModelStatusUnknown && return :Error
    error("Unknown HiGHS model status $(status)")
end
