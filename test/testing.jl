using Test, LinearAlgebra

include("../src/GMRES.jl")
include("../src/fea.jl")


@testset "hatfunction" begin
    mesh = meshproblem()

    ### Test at the middle of the elements
    y_test = [hatfunction(mesh[i], i, mesh) for i = 2:length(mesh)-1]
    @test y_test==ones(48)

    ### Test in the domain
    x_test = (mesh[2]+mesh[1])/2
    y_test = hatfunction(x_test, 2, mesh)
    y_gold = 0.5
    @test isapprox(y_test, y_gold)

    x_test = (mesh[2]+mesh[3])/2
    y_test = hatfunction(x_test, 2, mesh)
    @test isapprox(y_test, y_gold)

    ### Test the edges of the domain.
    #Test the lower edge 
    x_test = mesh[1]
    y_test = hatfunction(x_test, 2, mesh)
    y_gold = 0
    @test isapprox(y_test, y_gold)

    #Test that just inside the domain is infact inside the domain
    x_test = mesh[1]+eps()
    y_test = hatfunction(x_test, 2, mesh)
    y_gold = 0
    @test isapprox(y_test, y_gold, atol=1e-12)
    @test y_test>0

    #Test the upper edge
    x_test = mesh[3]
    y_test = hatfunction(x_test, 2, mesh)
    y_gold = 0
    @test isapprox(y_test, y_gold)

    #Test that just inside the domain is infact inside the domain
    x_test = mesh[3]-eps()
    y_test = hatfunction(x_test, 2, mesh)
    y_gold = 0
    @test isapprox(y_test, y_gold, atol=1e-12)
    @test y_test>0


    ### Test outside of the domain
    x_test = mesh[3]+eps()
    y_test = hatfunction(x_test, 2, mesh)
    y_gold = 0
    @test y_test==y_gold

    x_test = 0.75
    y_test = hatfunction(x_test, 2, mesh)
    y_gold = 0
    @test y_test==y_gold

    ### Test that you can't do basis functions on the first or last element. 
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

