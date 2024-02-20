using LinearAlgebra

"""
    mygmres(iters, b, x0, n, M, A)
Use the GMRES algorithm to solve a linear system of equations. 

**Inputs**
- iters::Int - The maximum number of iterations to take.
- b::Vector{Float} - The constants of the system. 
- x0::Vector{Float} - The initial guess of the solution.
- n::Int -  The dimension of the problem. 
- M::Matrix{Float} - The left-preconditioner matrix.
- A::Matrix{Float} - The coefficients of the system. 
"""
function mygmres(iters, b, x0, n, M, A; atol=1e-6)
    ### Preconditioning
    M_inv = inv(M)
    A = M_inv*A
    b = M_inv*b

    ### Initialize data storage. 
    v = zeros(length(b),iters+1)
    omega = zeros(length(b))
    h = zeros(iters+1,iters)

    ### Initial errors and Krylov vector. 
    r0 = b - A*x0
    beta = norm(r0)
    v[:,1] = r0/beta

    ### Arnoldi Iteration to get the Krylov vectors. 
    for j=1:iters
        omega .= A*v[:,j]

        for i=1:j
            h[i,j] = dot(omega,v[:,i])
            omega -= h[i,j]*v[:,i]
        end
        h[j+1,j] = norm(omega)

        #Avoid dividing by zero
        if h[j+1,j] == 0 #Todo: Should probably be a isapprox
            iters = j
            break
        end

        # Calculate the next Krylov vector
        v[:,j+1] .= omega/h[j+1,j]
    end

    ### Error vector? 
    e1 = zeros(iters+1)
    e1[1] = 1

    ### Calculate the result
    ym = h\(beta*e1)
    xm = x0 + v[:,1:end-1]*ym
    
    return xm
end