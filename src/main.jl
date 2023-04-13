"""
Run a simulation of the model.
"""
function main(;n=8000, H=0.08, R=1e-3, nstr=5, L=4e-2, eb=2.5e6, start=5e-3, tend=15e-9)
    T = Float64
    try
        plt.matplotlib.pyplot.style.use("granada")
    catch exc
        @warn "Unable to load matplotlib style.  This is fine: pictures will be uglier."
    end
    
    params = Parameters(;n, H, R, nstr, L, eb)
    
    @info("The following parameters are used in the simulation:\n\n" *
          sprint(io -> pretty_print(io, params, 1), context=:color => true))
    
    
    aux = AuxFields{T}(params)

    alpha = L^2 / (nstr * π * R^2)
    
    q0 = zeros(n)
    ne0 = zeros(n + 1)
    ztip0 = zeros(1)    
    
    ztip0[] = start
    ne0 .= @. (params.n0 * erfc((aux.zf - ztip0[]) / params.delta) / 2)
    
    u0 = ArrayPartition(q0, ne0, ztip0)

    prob = ODEProblem(derivs!, u0, (0.0, tend), (aux, params))

    @info "Solving the ODE."
    sol = solve(prob, Midpoint(), dtmax=1e-13, saveat=1e-9)
    @info "Done."
    
    return NamedTuple(Base.@locals)
end
