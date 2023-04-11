struct NaidisFit{T <: Number}
    c::T
    alpha::T
    e0::T
    c1::T
    alpha1::T
    e01::T
    v0::T
end

"Fit for positive streamer velocities. "
const fit_positive = NaidisFit(0.00172218,   # c
                               1.74504,      # alpha
                               6.10718e+06,  # e0
                               -0.000195532, # c1
                               0.648182,     # alpha1
                               6.10718e+06,  # e01
                               0.0)            # v0 

"Fit for negative streamer velocities. "
const fit_negative = NaidisFit(1.92409e-06,  # c
                               2.14334,      # alpha
                               4.3691e+06,   # e0
                               0.0279808,    # c1
                               1.0,          # alpha1
                               0.0,          # e01
                               44502.3)      # v0

"""
`naidis(fit, radius, ehead)`

Calculates the velocity of a streamer with given `radius` and peak field 
`ehead`
"""
function velocity(fit::NaidisFit, ehead, radius)
    max_velocity = 3e8
    ehead > fit.e0 || return 0.0
    
    v = (fit.c * (ehead - fit.e0)^fit.alpha * radius +
         fit.c1 * (ehead - fit.e01)^fit.alpha1 + fit.v0)

    min(v, max_velocity)
end


"""
`naidis(polarity, radius, ehead)`

Calculates the velocity of a streamer with given `radius` and peak field 
`ehead`
"""
function velocity(polarity::Integer, ehead, radius)
    tpl = (fit_negative, fit_positive)
    i = 1 + (polarity + 1) ÷ 2
    velocity(tpl[i], ehead, radius)
end
