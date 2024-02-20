


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
        if x<=mesh[i]
            return (x-mesh[i-1])/(mesh[i]-mesh[i-1])
        elseif x<mesh[i+1]
            return 1 - (x-mesh[i])/(mesh[i+1]-mesh[i])
        end
        return nothing
    else #Outside the domain of the basis function. 
        return 0
    end
end


"""
    meshproblem(; domain=(0,1), n_elems=50)

Create the vector across the domain. (The mesh is a vector across the domain). 
"""
function meshproblem(;domain=(0,1), n_elems=50)
    mesh = range(domain..., length=n_elems )
    return collect(mesh)
end