import Clp
import Cbc

function solve_IP(::Val{:cbc}, O::OptProblem, initial_solution = Int[],
                  use_warmstart = true; solver_options, solver_time, log_level)
    model = Cbc.Cbc_newModel()
    Cbc.Cbc_setParameter(model, "logLevel", string(max(0, log_level)))
    Cbc.Cbc_setParameter(model, "seconds", string(solver_time))
    Cbc.Cbc_setParameter(model, "allowableGap", "0")

    for (name, value) in solver_options
        Cbc.Cbc_setParameter(model, string(name), string(value))
    end

    A, lb, ub = add_cycle_constraints_to_formulation(O)
    cbc_loadproblem!(model, A, O.l, O.u, O.c, lb, ub)
    Cbc.Cbc_setObjSense(model, 1) # Minimize

    for i in 1:size(A, 2)
        Cbc.Cbc_setInteger(model, i - 1)
    end
    if !isempty(initial_solution) && use_warmstart
        Cbc.Cbc_setMIPStartI(model, length(initial_solution),
                             collect(Cint.(eachindex(initial_solution) .- 1)),
                             Float64.(initial_solution))
    end
    Cbc.Cbc_solve(model)

    attrs = Dict()
    attrs[:objbound] = Cbc.Cbc_getBestPossibleObjValue(model)
    attrs[:solver] = :ip
    solution = Solution(cbc_status(model), Cbc.Cbc_getObjValue(model),
                        copy(unsafe_wrap(Array,
                                         Cbc.Cbc_getColSolution(model),
                                         (size(A, 2),))),
                        attrs)
    Cbc.Cbc_deleteModel(model)
    return solution
end

function cbc_loadproblem!(model, A, l, u, c, lb, ub)
    Cbc.Cbc_loadProblem(model, size(A, 2), size(A, 1),
                        Cbc.CoinBigIndex.(A.colptr .- 1),
                        Int32.(A.rowval .- 1),
                        Float64.(A.nzval),
                        Float64.(l), Float64.(u), Float64.(c),
                        Float64.(lb), Float64.(ub))
end

function cbc_status(model)
    for (predicate, value) in ((Cbc.Cbc_isProvenOptimal, :Optimal),
                               (Cbc.Cbc_isProvenInfeasible, :Infeasible),
                               (Cbc.Cbc_isContinuousUnbounded, :Unbounded),
                               (Cbc.Cbc_isNodeLimitReached, :UserLimit),
                               (Cbc.Cbc_isSecondsLimitReached, :UserLimit),
                               (Cbc.Cbc_isSolutionLimitReached, :UserLimit),
                               (Cbc.Cbc_isAbandoned, :Error))
        predicate(model) != 0 && return value
    end
    return :InternalError
end
