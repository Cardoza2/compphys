using Plots

include("../src/GMRES.jl")
include("../src/fea.jl")

lambda = 1.

function u_ana(x)
    c1 = 1/(1-exp(1))
    c2 = -c1

    return c1*exp(x) + c2 + x
end

function print_mat(A)
    n, _ = size(A)
    for i in 1:n
        println(A[i,:])
    end
end


mesh = meshproblem(;n_nodes=50)


A = generate_A(lambda, mesh) 
b = generate_b(lambda, mesh)

# display(A)
# display(b)

n = length(b)

# u_weights = A\b
u_guess = zeros(n)
# u_guess = rand(n)
M = Diagonal(ones(n))

u_weights, err = mygmres(n, b, u_guess, n, M, A; tolerance=1e-8)

u(x) = evaluate_basis(x, mesh; weights=u_weights)

u_test = u.(mesh)
u_gold = u_ana.(mesh)


plt = plot(mesh, u_test, lab="FEA")
plot!(plt, mesh, u_gold, lab="Analytical")
display(plt)


# display(u_test)
# display(u_gold)

nothing