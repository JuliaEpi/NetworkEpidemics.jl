
"""
    meanfield(em, x0, kwargs)

Solve the mean field equations for the `AbstractEpidemicModel` `em`, from initial condition `x0`.

Keyword argument:
* `tmax=100.0`: Maximum time of the solution.
* `saveat=[]`: Times at which the solution should be computed. See the `DifferentialEquations` documentation on `solve` for more information.
"""
function meanfield(em::AbstractEpidemicModel, state_0; tmax=100.0, saveat=[])
    f! = meanfield_fun(em)
    u0 = init_state_mf(em, state_0)
    prob = ODEProblem(f!, u0, (0.0, tmax))
    sol = solve(prob, Tsit5(), saveat=saveat)
    finalize_meanfield(em, sol)
end

init_state_mf(em::AbstractEpidemicModel, x0) = float(x0)

finalize_meanfield(em::AbstractEpidemicModel, sol) = sol.t, sol.u, sol
