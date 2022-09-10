using DIVEnamics
using Test
using Random

@testset "DIVEnamics.jl" begin
    # Write your tests here.
end

points = 1000 * rand(MersenneTwister(42), 3, 1000)
point_attrs = randn(MersenneTwister(0), 1, 1000)
simplex_indices, simplex_nbs = delaunay(points)
@time void_cat = tetrahedra_to_voids(points, simplex_indices, false)
local_density = dtfe(points, simplex_indices, simplex_nbs)
println(void_cat[:,1])

#points = 1000 * rand(3, 1000000)
#simplex_indices = delaunay(points)
#@time void_cat = tetrahedra_to_voids(points, simplex_indices, false)
#mask = DIVEnamics.mask_in_box(void_cat, [0.,0.,0.], [1000., 1000., 1000.])

