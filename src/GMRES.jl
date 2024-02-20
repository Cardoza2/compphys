using LinearAlgebra

function mygmres(I, b, x0, n, M, A)
    M_inv = inv(M)
    A = M_inv*A
    b = M_inv*b

    v = zeros(length(b),I+1)
    ω = zeros(length(b))
    h = zeros(I+1,I)

    r0 = b - A*x0
    β = norm(r0)
    v[:,1] = r0/β

    for j=1:I
        ω .= A*v[:,j]
        for i=1:j
            h[i,j] = dot(ω,v[:,i])
            ω -= h[i,j]*v[:,j]
        end
        h[j+1,j] = norm(ω)
        if h[j+1,j] == 0
            I = j
            break
        end
        v[:,j+1] = ω/h[j+1,j]
    end
    e1 = zeros(length(b))
    e1[1] = 1
    ym = h\(β*e1)
    xm = x0 + Vm*ym
    return xm
end

A = rand(4,4)
b = rand(4)
x0 = rand(4)

first = A\b
second = mygmres(10,b,x0,4,I(4),A)
