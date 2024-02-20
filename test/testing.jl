using Test, LinearAlgebra

include("../src/GMRES.jl")


@testset "hatfunction" begin

end

@testset "GMRES" begin

    ### Example 1 in Fan reading
    A = [1.0 4 7; 2 9 7; 5 8 3]
    b = [1; 8; 2]

    x0 = [0,0,0]
    xstar =  mygmres(3, b, x0, 3, I(3), A)
    xtrue = A\b

    @test isapprox(xstar, xtrue)
    
    
    ### Compare against 
    n = 20
    A = rand(n,n)
    b = rand(n)
    x0 = zeros(n)
    iters = 20

    err =  max((mygmres(iters,b,x0,n,I(n),A) .- A\b)...)
    @test isapprox(err, 0, atol=1e-12)
    
end

