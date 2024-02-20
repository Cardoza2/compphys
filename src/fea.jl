
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




