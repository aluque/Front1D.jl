module Front1D

using Constants: co
using RecursiveArrayTools
using SpecialFunctions
using LoopVectorization
using Polyester
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using QuadGK
using PyPlot

const E_MOBILITY = 0.0372193
const IONIZATION_ALPHA = 433200.0
const IONIZATION_FIELD = 2e7
const ATTACHMENT_ALPHA = 2000.0
const ATTACHMENT_FIELD = 3e6

function main(;n=8000, H=0.08, R=.3e-3, nstr=5, L=4e-2, eb=2.5e6, start=2e-3)
    T = Float64
    plt.matplotlib.pyplot.style.use("granada")

    params = Parameters(;n, H, R, nstr, L, eb)

    
    aux = AuxFields{T}(params)

    @info "Pre-computed data finalized."
    alpha = L^2 / (nstr * π * R^2)
    
    q0 = zeros(n)
    ne0 = zeros(n + 1)
    ztip0 = zeros(1)    
    
    ztip0[] = start
    ne0 .= @. (params.n0 * erfc((aux.zf - ztip0[]) / params.delta) / 2)
    
    u0 = ArrayPartition(q0, ne0, ztip0)

    prob = ODEProblem(derivs!, u0, (0.0, 30e-9), (aux, params))
    @info "Solving the ODE."
    
    sol = solve(prob, Midpoint(), dtmax=1e-13, saveat=1e-9)

    @info "Done."
    
    return NamedTuple(Base.@locals)
end

@kwdef struct Parameters{T}
    "Number of cells"
    n::Int

    "Domain length"
    H::T
 
    "Cell length"
    dz::T = H / n
    
    "Streamer radius"
    R::T

    "Background field"
    eb::T
    
    "Field over-screening length"
    k::T = 1 / R
    
    "Number of streamers"
    nstr::Int

    "Cross-sectional area"
    L::T

    "Fraction of area covered by streamers"
    alpha::T = L^2 / (nstr * π * R^2)

    "Artificial length-scale for ionization"
    delta::T = R / 4

    "Constant velocity (to be improved)"
    v::T = 1.3e6

    "Electron density inside the streamer"
    neinside::T = 5e20
    
    "Constant n0 (to improve)"
    n0::T = neinside * nstr * π * R^2 / L^2
end


"""
A struct to store pre-allocated auxiliary fields.
"""
struct AuxFields{T, SB, F0, F1}
    zc::LinRange{T, Int}
    zf::LinRange{T, Int}
    
    jz::Vector{T}
    E::Vector{T}
    phi0::Vector{T}
    phi1::Vector{T}
    segbuf::SB
    
    A0::F0
    A1::F1
end

"""
Initialize an `AuxFields` instance with space for `n` cells and type `T`.
"""
function AuxFields{T}(params) where T
    (;n, dz, k, H) = params
    
    zc = LinRange{T}(dz / 2, H - dz / 2, n)
    zf = LinRange{T}(0, H, n + 1)

    jz = zeros(T, n + 1)
    E  = zeros(T, n + 1)
    phi0 = zeros(T, n)
    phi1 = zeros(T, n)

    segbuf = alloc_segbuf()
    A0 = ematrix(n, dz, 0)
    A1 = ematrix(n, dz, k)
    
    AuxFields{T, typeof(segbuf), typeof(A0), typeof(A1)}(zc, zf, jz, E, phi0, phi1, segbuf, A0, A1)
end

function ematrix(n, dz, k)
    (dz, k) = promote(dz, k)
    A = zeros(typeof(dz), n, n)
    
    for i in 1:n
        A[i, i] = -2 * co.epsilon_0 / dz^2 - k^2 * co.epsilon_0
        if i < n
            A[i, i + 1] = co.epsilon_0 / dz^2
        else
            # Dirichlet
            A[i, i] -= co.epsilon_0 / dz^2
        end
        
        if i > 1
            A[i, i - 1] = co.epsilon_0 / dz^2
        else
            # Dirichlet
            A[i, i] -= co.epsilon_0 / dz^2
        end
        
    end

    return factorize(A)
end


function poisson!(aux, params, q)
    (;n, dz, eb, alpha) = params
    
    (;phi0, phi1, E, A0, A1) = aux

    n = length(q)
    
    ldiv!(phi0, A0, q)
    ldiv!(phi1, A1, q)

    for i in 2:n
        E[i] = eb + (phi0[i] + (alpha - 1) * phi1[i] - phi0[i - 1] - (alpha - 1) * phi1[i - 1]) / dz
    end
    E[1] = eb + 2 * (phi0[1] + (alpha - 1) * phi1[1]) / dz
    E[n + 1] = eb - 2 * (phi0[n] + (alpha - 1) * phi1[n]) / dz    
end


function derivs!(du, u, (aux, params), t)
    (q, ne, ztip) = u.x
    (dq, dne, dztip) = du.x

    (;n, dz, delta, n0, v, nstr, R, L, eb) = params
    (;E, zc) = aux
    
    poisson!(aux, params, q)

    iztip = searchsortedlast(zc, ztip[])
    f = (ztip[] - zc[iztip]) / dz
    qtip = q[iztip] * (1 - f) + q[iztip + 1] * f
    Etip = qtip * L^2 / (2π * co.epsilon_0 * nstr * R) + eb
    ne1 = ionizationint(Etip) * nstr * π * R^2 / L^2

    #@show t Etip ne1
    
    @tturbo for i in 1:n
        dq[i] = co.elementary_charge * E_MOBILITY * (E[i] * ne[i] - E[i + 1] * ne[i + 1]) / dz
    end

    @tturbo for i in 1:n + 1
        z = dz * (i - 1)
        dne[i] = (effective_nu(E[i]) * ne[i] * erfc((z - ztip[1]) / delta) / 2
                  + ne1 * v * exp(-(z - ztip[1])^2 / delta^2) / (delta * sqrt(π)))
    end

    dztip[1] = v
end


impact_nu(eabs) = (E_MOBILITY * eabs * IONIZATION_ALPHA * exp(-IONIZATION_FIELD / eabs))
townsend(eabs) = IONIZATION_ALPHA * exp(-IONIZATION_FIELD / eabs)
attachment_nu(eabs) = (E_MOBILITY * eabs * ATTACHMENT_ALPHA * exp(-ATTACHMENT_FIELD / eabs))
effective_nu(eabs) = impact_nu(eabs) - attachment_nu(eabs)

function ionizationint(emax)
    return co.epsilon_0 * quadgk(townsend, 0.0, emax)[1] / co.elementary_charge
end


# ========================
# PLOT FUNCTIONS
# ========================
function plot_field(sol; kw...)
    (aux, params) = sol.prob.p
    (;H, n) = params
    (;zf) = aux
    
    for u in sol.u
        (q, ne, ztip) = u.x

        poisson!(aux, params, q)
        plt.plot(zf, aux.E; kw...)
    end        
end

function plot_q(sol; kw...)
    (aux, params) = sol.prob.p
    (;H, n) = params
    (;zc) = aux

    for u in sol.u
        (q, ne, ztip) = u.x

        plt.plot(zc, q; kw...)
    end        
end

function plot_ne(sol; kw...)
    (aux, params) = sol.prob.p
    (;H, n) = params
    (;zf) = aux

    for u in sol.u
        (q, ne, ztip) = u.x

        plt.plot(zf, ne; kw...)
    end        
end

end # module Front1D
