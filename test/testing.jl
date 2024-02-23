using Test, LinearAlgebra, Plots

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
    try
        y_test = hatfunction(x_test, 1, mesh)
    catch err_message
        @test err_message.msg=="hatfunction(): i should never equal 1."
    end

    try
        y_test = hatfunction(x_test, length(mesh), mesh)
    catch err_message
        n_mesh = length(mesh) #50
        err_gold = "hatfunction(): i should never equal $n_mesh. "
        @test err_message.msg==err_gold
    end


    ### Plot basis functions. 
    mesh = meshproblem(;n_nodes=7)
    xvec = append!(collect(range(0, 1, 100)), mesh)
    sort!(unique!(xvec))

    ymat = zeros(length(xvec), length(mesh)-2)

    testplt = plot(xaxis="Domain", yaxis="Basis functions", leg=:right)
    for i = 2:length(mesh)-1
        basis_idx = i-1
        ymat[:, basis_idx] = hatfunction.(xvec, Ref(i), Ref(mesh))
        
        plot!(testplt, xvec, ymat[:, basis_idx], lab="Basis $basis_idx")
    end
    vline!(testplt, mesh, lab=false, lw=2, alpha=0.5, linecolor=:black)
    # display(testplt)
    # savefig(testplt, "BasisFunctions.png")

end

@testset "hatfunc_deriv" begin
    mesh = meshproblem()

    # ### Test at the middle of the elements
    # y_test = [hatfunc_deriv(mesh[i], i, mesh) for i = 2:length(mesh)-1]
    # @test y_test==ones(48)

    # ### Test in the domain
    # x_test = (mesh[2]+mesh[1])/2
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # y_gold = 0.5
    # @test isapprox(y_test, y_gold)

    # x_test = (mesh[2]+mesh[3])/2
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # @test isapprox(y_test, y_gold)

    # ### Test the edges of the domain.
    # #Test the lower edge 
    # x_test = mesh[1]
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # y_gold = 0
    # @test isapprox(y_test, y_gold)

    # #Test that just inside the domain is infact inside the domain
    # x_test = mesh[1]+eps()
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # y_gold = 0
    # @test isapprox(y_test, y_gold, atol=1e-12)
    # @test y_test>0

    # #Test the upper edge
    # x_test = mesh[3]
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # y_gold = 0
    # @test isapprox(y_test, y_gold)

    # #Test that just inside the domain is infact inside the domain
    # x_test = mesh[3]-eps()
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # y_gold = 0
    # @test isapprox(y_test, y_gold, atol=1e-12)
    # @test y_test>0


    # ### Test outside of the domain
    # x_test = mesh[3]+eps()
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # y_gold = 0
    # @test y_test==y_gold

    # x_test = 0.75
    # y_test = hatfunc_deriv(x_test, 2, mesh)
    # y_gold = 0
    # @test y_test==y_gold

    # ### Test that you can't do basis functions on the first or last element. 
    # try
    #     y_test = hatfunc_deriv(x_test, 1, mesh)
    # catch err_message
    #     @test err_message.msg=="hatfunc_deriv(): i should never equal 1."
    # end

    # try
    #     y_test = hatfunc_deriv(x_test, length(mesh), mesh)
    # catch err_message
    #     n_mesh = length(mesh) #50
    #     err_gold = "hatfunc_deriv(): i should never equal $n_mesh. "
    #     @test err_message.msg==err_gold
    # end


    ### Plot basis functions. 
    mesh = meshproblem(;n_nodes=7)
    xvec = append!(collect(range(0, 1, 100)), mesh)
    sort!(unique!(xvec))

    ymat = zeros(length(xvec), length(mesh)-2)

    testplt = plot(xaxis="Domain", yaxis="Basis functions", leg=:right)
    for i = 2:length(mesh)-1
        basis_idx = i-1
        ymat[:, basis_idx] = hatfunc_deriv.(xvec, Ref(i), Ref(mesh))
        
        plot!(testplt, xvec, ymat[:, basis_idx], lab="Basis $basis_idx")
    end
    # vline!(testplt, mesh, lab=false, lw=2, alpha=0.5, linecolor=:black)
    # display(testplt)
    # savefig(testplt, "BasisFunctions.png")

end

@testset "evaluate_basis" begin
    ### test if it returns zeros at the edge of the domain
    mesh = meshproblem(;n_nodes=6)
    xvec = append!(collect(range(0, 1, 100)), mesh, [mesh[2]/2, (mesh[end-1]+mesh[end])/2, 0.5])
    sort!(unique!(xvec))

    yshape = evaluate_basis.(xvec, Ref(mesh))

    @test yshape[1]==0
    @test yshape[end]==0

    ### Test the sloped regions
    idx = findfirst(i->i==mesh[2]/2, xvec)
    @test yshape[idx]==0.5

    idx = findfirst(i->i==(mesh[end-1]+mesh[end])/2, xvec)
    @test yshape[idx]==0.5

    ### Test interior region
    idx = findfirst(i->i==0.5, xvec)
    @test yshape[idx]==1


    ### Plot shape function. 
    ymat = zeros(length(xvec), length(mesh)-2)
    

    testplt = plot(xaxis="Domain", yaxis="Basis functions", leg=:right)
    for i = 2:length(mesh)-1
        basis_idx = i-1
        ymat[:, basis_idx] = hatfunction.(xvec, Ref(i), Ref(mesh))
        
        plot!(testplt, xvec, ymat[:, basis_idx], lab="Basis $basis_idx")
    end
    plot!(testplt, xvec, yshape, lab="Shape Function", lw=1.5, linestyle=:dash)
    vline!(testplt, mesh, lab=false, lw=2, alpha=0.5, linecolor=:black)
    # display(testplt)
    # savefig(testplt, "ShapeFunction.png")
end

@testset "GMRES" begin

    ### Example 1 in Fan reading
    A = [1.0 4 7; 2 9 7; 5 8 3]
    b = [1.; 8; 2]

    x0 = [0.,0,0]
    xstar, err =  mygmres(3, b, x0, 3, I(3), A)
    xtrue = A\b

    @test isapprox(xstar, xtrue)
    
    
    ### Compare against 
    n = 20
    A = rand(n, n)
    b = rand(n)
    x0 = zeros(n)
    iters = 20

    y_test, err = mygmres(iters, b, x0, n, I(n), A; tolerance=1e-12)
    y_gold = A\b
    err =  maximum(y_test .- y_gold)
    @test isapprox(err, 0., atol=1e-12)

    # @time mygmres(iters, b, x0, n, I(n), A; tolerance=1e-12)
    
end

