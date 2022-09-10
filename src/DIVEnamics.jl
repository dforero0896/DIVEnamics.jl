module DIVEnamics
using TetGen
using GeometricalPredicates
using StaticArrays
using Zygote

export delaunay, tetrahedra_to_voids
export dtfe

GeometricalPredicates.Point(vec) = Point3D(Float64(vec[1]), Float64(vec[2]), Float64(vec[3]))
GeometricalPredicates.Primitive(matrix::SMatrix{3,4}) = Primitive([Point(matrix[:,j]) for j in 1:4]...)

function delaunay(points::Matrix)
    println("Starting Delaunay Triangulation with ", size(points, 2), " points.")
    @time begin
        input=TetGen.RawTetGenIO{Cdouble}(pointlist=points)
        tmp = tetrahedralize(input, "n")
        return tmp.tetrahedronlist, tmp.neighborlist
    end
end

GeometricalPredicates.volume(matrix::SMatrix{3,4}) = GeometricalPredicates.volume(GeometricalPredicates.Primitive(matrix))
circumcenter_x(matrix::SMatrix{3,4}) = circumcenter(GeometricalPredicates.Primitive(matrix))._x
circumcenter_y(matrix::SMatrix{3,4}) = circumcenter(GeometricalPredicates.Primitive(matrix))._y
circumcenter_z(matrix::SMatrix{3,4}) = circumcenter(GeometricalPredicates.Primitive(matrix))._z
circumradius(matrix::SMatrix{3,4}) = sqrt(circumradius2(GeometricalPredicates.Primitive(matrix)))

function mask_in_box(array::Matrix{<:AbstractFloat}, minima::Vector{<:AbstractFloat}, maxima::Vector{<:AbstractFloat})
    [all((array[1:3, i] .> minima) .& (array[1:3, i] .< maxima)) for i in 1:size(array,2)]
end


function tetrahedra_to_voids(points::Matrix{Float64}, tetrahedra_indices::Matrix{Int32}, derivatives::Bool)

    println("Found ", size(tetrahedra_indices, 2), " simplices")
    println("Building void catalog.")

    final_cat = zeros(5, size(tetrahedra_indices,2))
    if derivatives
        println("Allocating data for derivatives.")
        dx_dpos = Array{Float64, 3}(undef, (3, 4, size(tetrahedra_indices,2)))
        dy_dpos = Array{Float64, 3}(undef, (3, 4, size(tetrahedra_indices,2)))
        dz_dpos = Array{Float64, 3}(undef, (3, 4, size(tetrahedra_indices,2)))
        dvol_dpos = Array{Float64, 3}(undef, (3, 4, size(tetrahedra_indices,2)))
        dr_dpos = Array{Float64, 3}(undef, (3, 4, size(tetrahedra_indices,2)))
    end

    for i in 1:size(tetrahedra_indices,2)
        
        vertex_coordinates = SMatrix{3,4}(points[:,tetrahedra_indices[:,i]])
        simplex = GeometricalPredicates.Primitive(vertex_coordinates)
        center = circumcenter(simplex)
        radius = sqrt(circumradius2(simplex))
        vol = volume(simplex)
        final_cat[1,i] = center._x
        final_cat[2,i] = center._y
        final_cat[3,i] = center._z
        final_cat[4,i] = radius        
        final_cat[5,i] = vol
        if derivatives
            view(dx_dpos, :,:,i) .= gradient(circumcenter_x, vertex_coordinates)[1]
            view(dy_dpos, :,:,i) .= gradient(circumcenter_y, vertex_coordinates)[1]
            view(dz_dpos, :,:,i) .= gradient(circumcenter_z, vertex_coordinates)[1]
            view(dvol_dpos, :,:,i) .= gradient(GeometricalPredicates.volume, vertex_coordinates)[1]
            view(dr_dpos, :,:,i) .= gradient(circumradius, vertex_coordinates)[1]
        end
        
        
    end

    if derivatives
        return final_cat, dx_dpos, dy_dpos, dz_dpos, dvol_dpos, dr_dpos
    else
        return final_cat
    end

    
end

function dtfe(points::Matrix{Float64}, tetrahedra_indices::Matrix{Int32}, neighbor_list::Matrix{Int32})
    
    dtfe_arr = zeros(size(points,2))
    for i in 1:size(tetrahedra_indices,2)

        vertex_coordinates = SMatrix{3,4}(points[:,tetrahedra_indices[:,i]])
        vol = GeometricalPredicates.volume(vertex_coordinates)
        
        for j in 1:4
            dtfe_arr[tetrahedra_indices[j,i]] += vol
        end


    end

    return 4 ./ dtfe_arr

end




#Module end
end

