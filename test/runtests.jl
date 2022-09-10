using DIVEnamics
using Test

@testset "DIVEnamics.jl" begin
    # Write your tests here.
end

points = 1000 * rand(3, 1000)
simplex_indices = delaunay(points)
@time void_cat = tetrahedra_to_voids(points, simplex_indices, true)
println(void_cat[:,1])
#points = 1000 * rand(3, 1000000)
#simplex_indices = delaunay(points)
#@time void_cat = tetrahedra_to_voids(points, simplex_indices, false)
#mask = DIVEnamics.mask_in_box(void_cat, [0.,0.,0.], [1000., 1000., 1000.])

