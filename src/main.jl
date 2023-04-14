"""
Run simulations of the model and produce the plots used in the paper.
"""
function main(;results=[])

    try
        plt.matplotlib.pyplot.style.use("granada")
    catch exc
        @warn "Unable to load matplotlib style.  This is fine but pictures will be uglier."
    end
    
    number_of_streamers = [2, 5]
    @info "Running simulations with" number_of_streamers
    
    if isempty(results)
        # Passing results we avoid having to rerun the simulations.
        results = map(nstr -> simulate(Parameters(;nstr)), number_of_streamers)
    end
    
    paper_figure(number_of_streamers, results)
    

    # return everything in scope
    return NamedTuple(Base.@locals())
end


"""
Run a simulation of the model.
"""
function simulate(params)
    T = Float64
    
    printstyled("Running a simulation with the following parameters:\n\n", color=:light_green)
    pretty_print(params, 1)

    # Pre-allocate axiliary fields
    aux = AuxFields{T}(params)
    
    # Initial conditions
    q0 = zeros(params.n)
    ne0 = zeros(params.n + 1)
    ztip0 = zeros(1)        
    ztip0[] = params.start
    ne0 .= @. (params.nseed * erfc((aux.zf - ztip0[]) / params.delta) / 2)
    
    # Combine all variables in a single vector
    u0 = ArrayPartition(q0, ne0, ztip0)

    prob = ODEProblem(derivs!, u0, (0.0, params.tend), (aux, params))

    local sol
    with_logger(TerminalLogger()) do
        sol = solve(prob, Midpoint(), dtmax=1e-13, saveat=10e-9, progress=true)
    end
    
    @info "Done."
    
    return sol
end
