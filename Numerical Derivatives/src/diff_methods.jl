using FiniteDifferences

function richardson(f,x)
    dx = extrapolate_fdm(central_fdm(2, 1), f, x)[1]
    return dx
end
