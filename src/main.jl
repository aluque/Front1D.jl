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

    if isempty(results)
        # Passing results we avoid having to rerun the simulations.
        results = map(nstr -> simulate(Parameters(;nstr, R=1e-3)), number_of_streamers)
    end
    
    plotkw = Dict(2 => (c="k", lw=1.25, dashes=(5, 1)),
                  5 => (c="#ff5555", lw=1.25))
    
    plt.figure("Simplified model", figsize=(8, 6))

    # plot electric field
    plt.subplot(3, 1, 1)
    plt.subplots_adjust(top=0.95, hspace=0.1)
    for i in eachindex(number_of_streamers)
        nstr = number_of_streamers[i]
        sol = results[i]
        kw = plotkw[nstr]

        plot_field(sol; kw...)
    end
    plt.semilogy()
    plt.ylabel(L"$E_{int}$ (V/m)")
    noxticks()
    
    # plot charge
    plt.subplot(3, 1, 2)
    for i in eachindex(number_of_streamers)
        nstr = number_of_streamers[i]
        sol = results[i]
        kw = plotkw[nstr]

        plot_ne(sol; norm=true, kw...)
    end
    plt.ylabel(L"$\bar{n}_e$ (10$^{17}$ m$^{-3}$)")
    noxticks()
    
    # plot charge
    plt.subplot(3, 1, 3)
    for i in eachindex(number_of_streamers)
        nstr = number_of_streamers[i]
        sol = results[i]
        kw = plotkw[nstr]

        plot_q(sol; norm=true, kw...)
    end
    plt.ylabel(L"$\bar{q}$ (10$^{-3}$ Cm$^{-3}$)")
    plt.xlabel(L"$z$ (cm)")
    
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

    sol = solve(prob, Midpoint(), dtmax=1e-13, saveat=2e-9, progress=true)
    @info "Done."
    
    return sol
end
