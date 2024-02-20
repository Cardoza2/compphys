

function mygmres(i, b, x0, n, M, A)
    #Initialize r

    r0 = b - A*x0

    β = sqrt(sum(r0))

    v[1] = r0/β

    for j=1:I
        ω[j] = A*v[j]
        for i=1:j
            h[i,j] = sum(ω[j].*v[i])
            ω[j] -= h[i,j]*v[j]
        end
        h[j+1,j] = sqrt(sum(ω[j]))
        if h[j+1,j] == 0
            I = j
            break
        end
        v[j+1] = ω[j]/h[j+1,j]
    end
    
    ym = h\(β*e1)
    xm = x0 + Vm*ym

    return xm

end




function arnoldi(A, x0, m)

    return h, v
end