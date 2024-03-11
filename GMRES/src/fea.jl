using SparseArrays

"""
    meshproblem(; domain=(0,1), n_elems=50)

Create the vector across the domain. (The mesh is a vector across the domain). 


"""
function meshproblem(;domain=(0,1), n_nodes=50)
    mesh = range(domain..., length=n_nodes )
    return collect(mesh)
end

"""
    hatfunction(x, i, mesh)

Evaluate the ith hat function at x. 

i should range from 2, to M+1 (total mesh length is M+2)
"""
function hatfunction(x, i, mesh)

    if i==1
        error("hatfunction(): i should never equal 1.")
    elseif i==length(mesh)
        n_mesh = length(mesh)
        error("hatfunction(): i should never equal $n_mesh. ")
    end

    #Is x in the domain of the basis function? 
    if mesh[i-1]<=x<=mesh[i+1] 
        #Conditional to make a piecewise hat function. 
        if x<=mesh[i] #First half of the domain
            return (x-mesh[i-1])/(mesh[i]-mesh[i-1])
        elseif x<=mesh[i+1] #Second half of the domain
            return 1 - (x-mesh[i])/(mesh[i+1]-mesh[i])
        end

        @warn("hatfunction(): You've somehow entered a case that can't happen.")
        return 0 

    else #Outside the domain of the basis function. 
        return 0
    end
end

"""
    hatfunc_deriv(x, i, mesh)

Evaluate the derivative of the ith hat function (basis function). 

**Inputs**
- x::Float - The location to be evaluated
- i::Int - the ith basis function. 
- mesh::Vector{Float} - the domain to be evaluated. 
"""
function hatfunc_deriv(x, i, mesh)
    if i==1
        error("hatfunction(): i should never equal 1.")
    elseif i==length(mesh)
        n_mesh = length(mesh)
        error("hatfunction(): i should never equal $n_mesh. ")
    end

    #Is x in the domain of the basis function? 
    if mesh[i-1]<=x<=mesh[i+1] 
        #Conditional to make a piecewise hat function. 
        if x<=mesh[i] #First half of the domain
            return 1/(mesh[i]-mesh[i-1])
        elseif x<=mesh[i+1] #Second half of the domain
            return -1/(mesh[i+1]-mesh[i])
        end

        @warn("hatfunction(): You've somehow entered a case that can't happen.")
        return 0 

    else #Outside the domain of the basis function. 
        return 0
    end
end

"""

"""
function evaluate_basis(x, mesh; weights=ones(length(mesh)-2))

    if x>maximum(mesh)||x<minimum(mesh)
        error("evaluate_basis(): x outside of domain.")
    end

    if x==mesh[1]||x==mesh[end]
        return 0
    end

    ### Find the index which the point is located
    node_idx = findfirst(item -> item>x, mesh)
    n = length(mesh)

    if node_idx==2 # bottom basis
        return weights[1]*hatfunction(x, 2, mesh)

    elseif node_idx==n # Top basis
        return weights[end]*hatfunction(x, node_idx-1, mesh)

    else #x is in a region that has two basis functions. 
        y1 = weights[node_idx-2]*hatfunction(x, node_idx-1, mesh) #lower 
        y2 = weights[node_idx-1]*hatfunction(x, node_idx, mesh) #upper

        return y1+y2
    end
end


function A1(i, j, mesh)
    if i+1==j
        top = -(mesh[j]-mesh[i])
        bot = (mesh[i+1]-mesh[i])*(mesh[j]-mesh[j-1])
        return top/bot

    elseif i == j+1
        top = -(mesh[i]-mesh[j])
        bot = (mesh[i]-mesh[i-1])*(mesh[j+1]-mesh[j])
        return top/bot

    elseif i == j
        t1 = 1/(mesh[i]-mesh[i-1])
        t2 = 1/(mesh[i+1]-mesh[i])
        return t1 + t2

    else
        return 0
    end
end

function A2(i, j, mesh)
    if i == j+1
        top = (mesh[i]^2 - mesh[i-1]^2)/2 - mesh[j]*(mesh[i] - mesh[i-1])
        bot = (mesh[j+1]-mesh[j])*(mesh[i]-mesh[i-1])
        return 1 - top/bot

    elseif i+1 == j
        top = mesh[j-1]*(mesh[i+1]-mesh[i]) - (mesh[i+1]^2 - mesh[i]^2)/2
        bot = (mesh[j]-mesh[j-1])*(mesh[i+1]-mesh[i])
        return top/bot

    elseif i==j
        top1 = (mesh[i]^2 + mesh[i-1]^2)/2 - mesh[i]*mesh[i-1]
        bot1 = (mesh[i]-mesh[i-1])^2

        top2 = (mesh[i+1]^2 + mesh[i]^2)/2 - mesh[i]*mesh[i+1]
        bot2 = (mesh[i+1] - mesh[i])^2

        return top1/bot1 - 1 + top2/bot2
    else
        return 0
    end
end


function generate_A(lambda, mesh)
    n = length(mesh)-2

    A = spzeros(n, n)
    ### Find A
    for i = 1:n
        for j = 1:n
            A[i, j] = A1(j+1, i+1, mesh)
            A[i, j] += lambda*A2(j+1, i+1, mesh)
        end
    end

    return A
end

function bval(i, mesh)
    top1 = (mesh[i]^2 - mesh[i-1]^2)/2 - mesh[i-1]*(mesh[i]-mesh[i-1])
    bot1 = (mesh[i]-mesh[i-1])

    t2 = mesh[i+1]-mesh[i]

    top3 = (mesh[i+1]^2 - mesh[i]^2)/2 - mesh[i]*(mesh[i+1]-mesh[i])
    bot3 = (mesh[i+1]-mesh[i]) 

    return top1/bot1 + t2 - top3/bot3
end

function generate_b(lambda, mesh)
    n = length(mesh)-2

    b = zeros(n) #Doesn't need to be sparse since b is dense.
    for i = 1:n
        b[i] = bval(i+1, mesh)
    end

    return b
end