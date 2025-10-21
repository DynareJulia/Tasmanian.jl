using Plots
using Random
using Tasmanian
using Tasmanian.RecipesBase
using Test
using VisualRegressionTests

@testset verbose=true "Testing Tasmanian" begin
    @testset verbose=true "Basic testing" include("test_tasmanian.jl")
    @testset verbose=true "Testing examples" include("test_examples.jl")
    @testset verbose=true "Aqua test" include("Aqua.jl")
    @testset verbose=true "JET test" include("JET.jl")
end

nothing
