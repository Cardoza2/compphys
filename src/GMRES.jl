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
function mygmres(iters, b, x0, n, M, A; tolerance=1e-6, verbose::Bool=false)
    ### Preconditioning
    M_temp = zeros(size(M))
    M_temp .= M
    M_inv = inv(M_temp)
    A = sparse(M_inv*A)
    b = sparse(M_inv*b)

    x_sol = deepcopy(x0)

    ### Initialize data storage.
    n = length(b)
    v = zeros(n,iters+1) #Krylov vectors
    omega = zeros(n)
    h = zeros(iters+1,iters) #Hessenberg matrix
    r0 = zeros(n) #Residual vector
    ym = zeros(iters+1)

    e1 = zeros(iters+1) #Error vector?
    e1[1] = 1

    for i = 2:iters #todo: skipping some might be more efficient
        v_i = view(v, :, 1:i+1)
        h_i = view(h, 1:i+1, 1:i)
        e_i = view(e1, 1:i+1)
        ym_i = view(ym, 1:i)

        Arnoldi!(v_i, h_i, omega, r0, x0, A, b, i, i-1)

        ### Calculate the result
        beta = norm(r0)
        ym_i .= h_i\(beta*e_i)
        x_sol .= x0 + v_i[:,1:end-1]*ym_i

        err = norm(A*x_sol - b, Inf)

        if err<=tolerance
            return x_sol
        end
    end

    if verbose
        @warn("mygres(): Maximum number of iterations hit without converging.")
    end
    return x_sol
end

"""
    mygmres!(v, h, omega, ...)

The Arnoldi iteration algorithm to generate the Krylov vectors. Note it generates the vectors in place.

**Inputs**
"""
function Arnoldi!(v, h, omega, r0, x0, A, b, iters, current_iter)

    ### Initial errors and Krylov vector.
    r0 .= b - A*x0
    beta = norm(r0)
    v[:,1] = r0/beta

    ### Arnoldi Iteration to get the Krylov vectors.
    for j=current_iter:iters
        omega .= A*v[:,j]

        for i=1:j
            h[i,j] = dot(omega,v[:,i])
            omega -= h[i,j]*v[:,i]
        end
        h[j+1,j] = norm(omega)

        #Avoid dividing by zero
        if isapprox(h[j+1,j], 0)
            iters = j
            break
        end

        # Calculate the next Krylov vector
        v[:,j+1] .= omega/h[j+1,j]
    end
end
