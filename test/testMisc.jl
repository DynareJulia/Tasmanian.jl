using Tasmanian, Tasmanian_jll
using Test
using Plots

"""
    Check the consistency of derivatives in terms of dimensions and ordering.
"""
function checkDerivatives()
    func1(x) = 2.0 * x[1] .* x[1] + x[2] .* x[2] / 2.0
    grid1 = makeGlobalGrid(
        dimension = 2, outputs = 1, depth = 4, type = "iptotal", rule = "gauss-legendre")
    points1 = getNeededPoints(grid1)
    values1 = mapslices(func1, points1, dims = 1)
    loadNeededValues!(grid1, values1)
    grad1 = differentiate(grid1, [3.0, 4.0])
    @test size(grad1) == (2, 1)
    @test grad1 ≈ [12.0, 4.0]

    function func2(x)
        [2.0 * x[1] .* x[1] + x[2] .* x[2] / 2.0 + x[3] .* x[3], x[1] .* x[2] .* x[3]]
    end
    grid2 = makeGlobalGrid(
        dimension = 3, outputs = 2, depth = 4, type = "iptotal", rule = "gauss-legendre")
    points2 = getNeededPoints(grid2)
    values2 = mapslices(func2, points2, dims = 1)
    loadNeededValues!(grid2, values2)
    grad2 = differentiate(grid2, [1.0, 2.0, 3.0])
    @test size(grad2) == (3, 2)
    @test grad2 ≈ [4.0 2.0 6.0; 6.0 3.0 2.0]'
end

"""
    Miscelaneous tests that don't quite fit in other categories.
"""

@testset verbose=true "Testing plotting and other misc" begin
    @testset "checkDerivatives" checkDerivatives()
end
