using Test

@testset verbose = true "Testing Tasmanian" begin
    @testset verbose = true "Testing Core I/O" include("testBasicIO.jl")
    @testset verbose = true "Testing accelerated evaluate consistency" include("testAcceleration.jl")
    @testset verbose = true "Testing error handling" include("testExceptions.jl")
    @testset verbose = true "Testing core make/update grid" include("testMakeUpdate.jl")
    @testset verbose = true "Testing core refine grid" include("testRefinement.jl")
    @testset verbose = true "Testing core learning from random samples" include("testUnstructuredData.jl")
    @testset verbose = true "Testing miscelaneous" include("testMisc.jl")
end

nothing
