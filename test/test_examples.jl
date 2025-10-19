@testset "Example 1" begin
    include("../examples/example1.jl")
    @plottest plot(p) "refplots/example1.png"
end
@testset "Example 2" begin
    include("../examples/example2.jl")
    @plottest plot(p) "refplots/example2.png"
end
@testset "Example 3" begin
    include("../examples/example3.jl")
    @plottest plot(p) "refplots/example3.png"
end
@testset "example_sparse_grids_01" include("../examples/example_sparse_grids_01.jl")
@testset "example_sparse_grids_02" include("../examples/example_sparse_grids_02.jl")
@testset "example_sparse_grids_03" include("../examples/example_sparse_grids_03.jl")
@testset "example_sparse_grids_04" include("../examples/example_sparse_grids_04.jl")
@testset "example_sparse_grids_05" include("../examples/example_sparse_grids_05.jl")
