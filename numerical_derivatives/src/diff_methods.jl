using FiniteDifferences, FastChebInterp

function richardson(f,x)
    dx = extrapolate_fdm(central_fdm(2, 1), f, x)[1]
    return dx
end

function chebyshev(f, x; domain=(-1, 1), order=200)
    lb = domain[1]
    ub = domain[2]
    xpnts = chebpoints(order, lb, ub) 
    c = chebinterp(f.(xpnts), lb, ub)

    if length(x)==1
     _, gapprox = chebgradient(c, x)
    else
        gapprox = zero(x)

        for i in eachindex(x)
            _, g = chebgradient(c, x[i])
            gapprox[i] = g
        end
    end

    return gapprox
end

function forward_AD(f, x; optargs...)
    return ForwardDiff.derivative(f,x; optargs...)
end

function ddx(f, x, method; optargs...)
    if method in ("chebyshev", "Chebyshev", "c", "C")
        return chebyshev(f,x; optargs...)
    elseif method in ("richardson", "Richardson", "r", "R")
        return richardson(f, x; optargs...)
    elseif method in ("forward_AD", "forwarddiff", "AD", "ad", "a", "A")
        return forward_AD(f, x; optargs...)
    else
        error("differentiation method $method not recognized! Choose from ('chebyshev', 'richardson', 'AD').")
    end
end
