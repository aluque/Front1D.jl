# ========================
# PLOT FUNCTIONS
# ========================
function plot_field(sol; norm=false, kw...)
    (aux, params) = sol.prob.p
    (;H, n, nstr) = params
    (;zf) = aux

    if norm
        @warn "norm keyword is ignored"
    end
    
    for u in sol.u
        (q, ne, ztip) = u.x

        poisson!(aux, params, q)
        plt.plot(zf, aux.E; kw...)
    end        
end

function plot_q(sol; norm=false, kw...)
    (aux, params) = sol.prob.p
    (;H, n, nstr) = params
    (;zc) = aux

    for u in sol.u
        (q, ne, ztip) = u.x
        f = norm ? 1 / nstr : oneunit(1 / nstr)

        plt.plot(zc, f * q; kw...)
    end        
end

function plot_ne(sol; norm=false, kw...)
    (aux, params) = sol.prob.p
    (;H, n, nstr) = params
    (;zf) = aux

    for u in sol.u
        (q, ne, ztip) = u.x
        f = norm ? 1 / nstr : oneunit(1 / nstr)

        plt.plot(zf, f * ne; kw...)
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
    plt.plot(sol.t, ztip; kw...)
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
