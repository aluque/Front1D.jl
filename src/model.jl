const E_MOBILITY = 0.0372193
const IONIZATION_ALPHA = 433200.0
const IONIZATION_FIELD = 2e7
const ATTACHMENT_ALPHA = 2000.0
const ATTACHMENT_FIELD = 3e6


"""
A struct to contain the model parameters.

The default values of the parameters are those used in the paper.
"""
@kwdef struct Parameters{T}
    "Number of cells"
    n::Int = 8000

    "Domain length (m)"
    H::T = 8.192e-2
 
    "Side length (m)"
    L::T = 4.096e-2

    "Cell length (m)"
    dz::T = H / n
    
    "Streamer radius (m)"
    R::T = 1e-3

    "Background field (V/m)"
    eb::T = 2.5e6
    
    "Number of streamers"
    nstr::Int

    "Starting location of the front (ztip) (m)"
    start::T = 5e-3
    
    "Field over-screening decay rate (1/length) (1/m)"
    k::T = 1 / R
    
    "1 / fraction of area covered by streamers"
    alpha::T = L^2 / (nstr * π * R^2)

    "Artificial length-scale for ionization (m)"
    delta::T = R

    "Constant velocity (if selected) (m/s)"
    v::T = 1.3e6

    "Electron density inside the streamer (1/m3)"
    neinside::T = 5e21
    
    "Constant tip electron density (if selected) (1/m3)"
    n0::T = 0.0

    "Linear increase of the tip electron density (1/m4)"
    nslope::T = 5e18 * nstr

    "Electron density in the initial (seed) front (1/m3)"
    nseed::T = 1e16 * nstr
    
    "Set a constant electron density at the tip"
    fix_netip::Bool = true

    "Set a constant velocity"
    fix_v::Bool = false

    "Final simulation time (s)"
    tend::T = 40e-9
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


"""
Compute a matrix factorization to solve the Poisson / Helmholtz problem.
"""
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

"""
Calculate the electric field due to charge q and store it in aux.E.
"""
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


"""
Compute the derivatives of all variables in the model stored in `u` and store them in `du`.
Requires `aux` and `params` passed as a tuple in the third argument.  Time is received for 
compatibility the the DiffEq API but ignored.
"""
function derivs!(du, u, (aux, params), t)
    (q, ne, ztip) = u.x
    (dq, dne, dztip) = du.x

    (;n, dz, delta, n0, nslope, v, nstr, R, L, eb, fix_netip, fix_v) = params
    (;E, zc, zf) = aux
    
    poisson!(aux, params, q)

    izctip = searchsortedlast(zc, ztip[])
    f = (ztip[] - zc[izctip]) / dz
    qtip = q[izctip] * (1 - f) + q[izctip + 1] * f


    # Naive spherical charge.
    # Etip = qtip * L^2 / (2π * co.epsilon_0 * nstr * R) + eb
    Etip = maximum(E)

    # A variable electron density.
    ne1 = ionizationint(Etip) * nstr * π * R^2 / L^2
    netip = fix_netip ? (n0 + nslope * ztip[]) : ne1

    # A variable velocity
    v1 = fix_v ? v : velocity(+1, Etip, R)

    @tturbo for i in 1:n
        dq[i] = co.elementary_charge * E_MOBILITY * (E[i] * ne[i] - E[i + 1] * ne[i + 1]) / dz
    end

    @tturbo for i in 1:n + 1
        z = dz * (i - 1)
        
        dne[i] = (effective_nu(E[i]) * ne[i] * erfc((z - ztip[1]) / delta) / 2
                  + netip * v1 * exp(-(z - ztip[1])^2 / delta^2) / (delta * sqrt(π)))        
    end

    dztip[1] = v1
end


impact_nu(eabs) = (E_MOBILITY * eabs * IONIZATION_ALPHA * exp(-IONIZATION_FIELD / eabs))
townsend(eabs) = IONIZATION_ALPHA * exp(-IONIZATION_FIELD / eabs)
attachment_nu(eabs) = (E_MOBILITY * eabs * ATTACHMENT_ALPHA * exp(-ATTACHMENT_FIELD / eabs))
effective_nu(eabs) = impact_nu(eabs) - attachment_nu(eabs)

Base.@assume_effects :effect_free @inline function ionizationint(emax)
    return co.epsilon_0 * quadgk(townsend, 0.0, emax)[1] / co.elementary_charge
end
