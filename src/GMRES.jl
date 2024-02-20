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
            ω -= h[i,j]*v[:,i]
        end
        h[j+1,j] = norm(ω)
        if h[j+1,j] == 0
            I = j
            break
        end
        v[:,j+1] .= ω/h[j+1,j]
    end
    e1 = zeros(I+1)
    e1[1] = 1
    ym = h\(β*e1)
    xm = x0 + v[:,1:end-1]*ym
    return xm
end

function run_test()

    A = [1.0 4 7; 2 9 7; 5 8 3]
    b = [1; 8; 2]

    x0 = [0,0,0]
    return mygmres(3, b, x0, 3, I(3), A), A\b

end

function run_test_big(;n=20,iters=10)
    
    A = rand(n,n)
    b = rand(n)
    x0 = zeros(n)

    return max((mygmres(iters,b,x0,n,I(n),A) .- A\b)...)

end