function paper_figure(number_of_streamers, results)
    # plotkw = Dict(2 => (c="k", lw=1.25, dashes=(5, 1)),
    #               5 => (c="#ff5555", lw=1.25))
    
    plt.figure("Simplified model", figsize=(8, 8))
    plt.subplots_adjust(top=0.95, hspace=0.1, wspace=0.06, bottom=0.15)
   
    cm = plt.get_cmap("plasma_r")

    # plot electric field
    plt.subplot(3, 2, 1)
    plot_field(results[1]; cm)
    plt.title("N=2 streamers")
    plt.semilogy()
    plt.ylim([5e5, 2e7])
    plt.ylabel(L"$E_{int}$ (V/m)")
    noxticks()

    plt.subplot(3, 2, 2)
    plot_field(results[2]; cm)
    plt.title("N=5 streamers")
    plt.semilogy()
    plt.ylim([5e5, 2e7])
    noxticks()
    noyticks()
    
    # plot ne
    plt.subplot(3, 2, 3)
    plot_ne(results[1]; cm, norm=true)
    plt.ylabel(L"$\bar{n}_e$ (10$^{17}$ m$^{-3}$)")
    plt.ylim([-.5, 3])
    noxticks()

    plt.subplot(3, 2, 4)
    plot_ne(results[2]; cm, norm=true)
    plt.ylim([-.5, 3])
    noxticks()
    noyticks()
    
    # plot charge
    plt.subplot(3, 2, 5)
    plot_q(results[1]; cm, norm=true)
    plt.ylabel(L"$\bar{q}$ (10$^{-4}$ Cm$^{-3}$)")
    plt.ylim([-2, 9])
    plt.xlabel(L"$z$ (cm)")

    plt.subplot(3, 2, 6)
    lines = plot_q(results[2]; cm, norm=true)
    plt.ylim([-2, 9])
    noyticks()
    plt.xlabel(L"$z$ (cm)")

    times_ns = round.(Int, results[1].t ./ 1e-9)
    labels = ["$tns ns" for tns in times_ns]
    
    plt.legend(lines, labels, loc=(-0.85, -0.5),
               bbox_transform=plt.gcf().transFigure, ncol=length(times_ns))

    plt.savefig("output/simplemodel.pdf")
end



# ========================
# PLOT FUNCTIONS
# ========================
function plot_field(sol; norm=false, cm=nothing, kw...)
    (aux, params) = sol.prob.p
    (;H, n, nstr, tend) = params
    (;zf) = aux

    if norm
        @warn "norm keyword is ignored"
    end
    
    for (i, u) in enumerate(sol.u)
        (q, ne, ztip) = u.x

        poisson!(aux, params, q)
        kw1 = isnothing(cm) ? kw : (color = cm(sol.t[i] / tend)[1:3], kw...)
        plt.plot(zf ./ co.centi, aux.E; kw1...)
    end        
end

function plot_q(sol; norm=false, cm=nothing, kw...)
    (aux, params) = sol.prob.p
    (;H, n, nstr, tend) = params
    (;zc) = aux

    lines = Any[]
    
    for (i, u) in enumerate(sol.u)
        (q, ne, ztip) = u.x
        f = norm ? 1 / nstr : oneunit(1 / nstr)

        kw1 = isnothing(cm) ? kw : (color = cm(sol.t[i] / tend)[1:3], kw...)
        line = plt.plot(zc ./ co.centi, f * q / 1e-4; kw1...)
        push!(lines, line[1])
    end

    return lines
end

function plot_ne(sol; norm=false, cm=nothing, kw...)
    (aux, params) = sol.prob.p
    (;H, n, nstr, tend) = params
    (;zf) = aux

    for (i, u) in enumerate(sol.u)
        (q, ne, ztip) = u.x
        f = norm ? 1 / nstr : oneunit(1 / nstr)

        kw1 = isnothing(cm) ? kw : (color = cm(sol.t[i] / tend)[1:3], kw...)
        plt.plot(zf ./ co.centi, f * ne / 1e17; kw1...)
    end        
end

function plot_ztip(sol; kw...)
    (aux, params) = sol.prob.p
    (;H, n) = params
    (;zf) = aux

    ztip = map(sol.u) do u
        (q, ne, ztip) = u.x

        poisson!(aux, params, q)
        imax = argmax(aux.E)
        return zf[imax]
    end
    plt.plot(sol.t, ztip ./ co.centi; kw...)
end

function plot_v(sol; kw...)
    (aux, params) = sol.prob.p
    (;H, n) = params
    (;zf) = aux

    ztip = map(sol.u) do u
        (q, ne, ztip) = u.x

        poisson!(aux, params, q)
        imax = argmax(aux.E)
        return zf[imax]
    end
    t1 = @. (0.5 * (sol.t[begin:end - 1] + sol.t[begin + 1:end]))
    plt.plot(t1, diff(ztip) ./ diff(sol.t); kw...)
end

""" Remove xticks from the plot. """
function noxticks(ax=plt.gca())
    loc = ax.get_xticks()
    ax.set_xticklabels(["" for l in loc])
end

""" Remove yticks from the plot. """
function noyticks(ax=plt.gca())
    loc = ax.get_yticks()
    ax.set_yticklabels(["" for l in loc])
end
