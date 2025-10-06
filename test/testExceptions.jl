using Tasmanian
using Test

include("testCommon.jl")

function getSparseGridTests()
    @testset verbose = true "makeGlobalGrid" begin
        @test_throws "dimension" makeGlobalGrid(dimension = -1, outputs = 1,  depth = 4, type = "level", rule = "clenshaw-curtis")
        @test_throws "outputs" makeGlobalGrid(dimension = 2, outputs = -1,  depth = 4, type = "level", rule = "clenshaw-curtis")
        @test_throws "depth" makeGlobalGrid(dimension = 2, outputs = 1,  depth = -4, type = "level", rule = "clenshaw-curtis")
        @test_throws "type" makeGlobalGrid(dimension = 2, outputs = 1,  depth = 4, type = "wrong", rule = "clenshaw-curtis")
        @test_throws "rule" makeGlobalGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "clenshaw-wrong")
        @test_throws "weights" makeGlobalGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "clenshaw-curtis", anisotropic_weights=[1,2,3])
        @test_throws "limits" makeGlobalGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "clenshaw-curtis", anisotropic_weights=[1,2], level_limits = [1, 2, 3])
        @test makeGlobalGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "clenshaw-curtis", anisotropic_weights=[1,2], level_limits = [1, 2]) isa TasmanianSG
        @test makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "chebyshev") isa TasmanianSG      
    end
    @testset verbose = true "makeSequenceGrid" begin
        @test_throws "dimension" makeSequenceGrid(dimension = -1, outputs = 1,  depth = 4, type = "level", rule = "leja")
        @test_throws "outputs" makeSequenceGrid(dimension = 2, outputs = -1,  depth = 4, type = "level", rule = "leja")
        @test_throws "depth" makeSequenceGrid(dimension = 2, outputs = 1,  depth = -4, type = "level", rule = "leja")
        @test_throws "type" makeSequenceGrid(dimension = 2, outputs = 1,  depth = 4, type = "wrong", rule = "leja")
        @test_throws "rule" makeSequenceGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "weja")
        @test_throws "weights" makeSequenceGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "leja", anisotropic_weights=[1,2,3])
        @test_throws "limits" makeSequenceGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "leja", anisotropic_weights=[1,2], level_limits = [1, 2, 3])
        @test makeSequenceGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", rule = "leja", anisotropic_weights=[1,2], level_limits = [1, 2]) isa TasmanianSG
        @test makeSequenceGrid(dimension = 2, outputs = 2, depth = 2, type = "level", rule = "rleja") isa TasmanianSG
    end
    @testset verbose = true "makeLocalPolynomialGrid" begin
        @test_throws "dimension" makeLocalPolynomialGrid(dimension = -1, outputs = 1,  depth = 4, order = 2, rule = "localp")
        @test_throws "outputs" makeLocalPolynomialGrid(dimension = 2, outputs = -1,  depth = 4, order = 2, rule = "localp")
        @test_throws "depth" makeLocalPolynomialGrid(dimension = 2, outputs = 1,  depth = -4, order = 2, rule = "localp")
        @test_throws "order" makeLocalPolynomialGrid(dimension = 2, outputs = 1,  depth = 4, order = -2, rule = "localp")
        @test_throws "rule" makeLocalPolynomialGrid(dimension = 2, outputs = 1,  depth = 4,  order = 2, rule = "lowrong")
        @test_throws "limits" makeLocalPolynomialGrid(dimension = 2, outputs = 1,  depth = 4, order = 2, rule = "localp", level_limits = [1, 2, 3])
        @test makeLocalPolynomialGrid(dimension = 2, outputs = 1,  depth = 4,  rule = "localp", order = 2, level_limits = [1, 2]) isa TasmanianSG
    end
    @testset verbose = true "makeWaveletGrid" begin
        @test_throws "dimension" makeWaveletGrid(dimension = -1, outputs = 1,  depth = 4, order = 1)
        @test_throws "outputs" makeWaveletGrid(dimension = 2, outputs = -1,  depth = 4, order = 1)
        @test_throws "depth" makeWaveletGrid(dimension = 2, outputs = 1,  depth = -4, order = 3)
        @test_throws "order" makeWaveletGrid(dimension = 2, outputs = 1,  depth = 4, order = 2)
        @test_throws "limits" makeWaveletGrid(dimension = 2, outputs = 1,  depth = 4, order = 1, level_limits = [1, 2, 3])
        @test makeWaveletGrid(dimension = 2, outputs = 1,  depth = 4,  order = 1, level_limits = [1, 2]) isa TasmanianSG
    end
    @testset verbose = true "makeFourierGrid" begin
        @test_throws "dimension" makeFourierGrid(dimension = -1, outputs = 1,  depth = 4, type = "level")
        @test_throws "outputs" makeFourierGrid(dimension = 2, outputs = -1,  depth = 4, type = "level")
        @test_throws "depth" makeFourierGrid(dimension = 2, outputs = 1,  depth = -4, type = "level")
        @test_throws "type" makeFourierGrid(dimension = 2, outputs = 1,  depth = 4, type = "wrong")
        @test_throws "limits" makeFourierGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", level_limits = [1, 2, 3])
        @test_throws "weights" makeFourierGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", anisotropic_weights=[1, 2, 3], level_limits = [1, 2]) isa TasmanianSG
        @test makeFourierGrid(dimension = 2, outputs = 1,  depth = 4, type = "level", anisotropic_weights=[1, 2], level_limits = [1, 2]) isa TasmanianSG
    end
    @testset "loadNeededValues!" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 2, depth = 2, type = "level", rule = "rleja")
        @test loadNeededValues!(grid, zeros(2, 6)) isa Nothing
    end
    @testset "updateGlobalGrid!" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "global" updateGlobalGrid!(grid, depth = 1, type = "iptotal")
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "chebyshev")
        @test_throws "depth" updateGlobalGrid!(grid, depth = -1, type = "iptotal")
        @test_throws "type" updateGlobalGrid!(grid, depth = 4, type = "wrong")
        @test_throws "weights" updateGlobalGrid!(grid, depth = 4, type = "iptotal", anisotropic_weights=[1,2,3])
        @test_throws "limits" updateGlobalGrid!(grid, depth =4, type = "iptotal", level_limits = [1,2,3])
        @test updateGlobalGrid!(grid, depth = 4, type = "iptotal", level_limits = [1,2]) isa Nothing
    end
    @testset "updateSequenceGrid!" begin
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "sequence" updateSequenceGrid!(grid, depth = 4, type = "iptotal")
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "depth" updateSequenceGrid!(grid, depth = -1,type = "iptotal")
        @test_throws "type" updateSequenceGrid!(grid, depth = 4, type = "wrong")
        @test_throws "weights" updateSequenceGrid!(grid, depth = 4, type = "iptotal", anisotropic_weights = [1,2,3])
        @test_throws "limits" updateSequenceGrid!(grid, depth = 4, type = "iptotal", level_limits = [1,2,3])
        @test updateSequenceGrid!(grid, depth = 4, type = "iptotal",level_limits = [2,3]) isa Nothing
    end
    @testset "updateFourier!" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "Fourier" updateFourierGrid!(grid, depth = 3, type = "level")
        grid = makeFourierGrid(dimension = 2, outputs = 1, depth = 2, type = "level")
        @test_throws "depth" updateFourierGrid!(grid, depth = -3, type = "level")
        @test_throws "type" updateFourierGrid!(grid, depth = 3, type ="wrong")
        @test_throws "weight" updateFourierGrid!(grid, depth = 3, type = "iptotal", anisotropic_weights=[1])
        @test_throws "limits" updateFourierGrid!(grid, depth = 3, type = "iptotal", level_limits = [4])
        @test updateFourierGrid!(grid, depth = 3, type = "iptotal", anisotropic_weights=[4, 3], level_limits=[4, 3]) isa Nothing
    end
    @testset "getInterpolationWeights" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "length(x)" getInterpolationWeights(grid, [1,2,3])
        @test_throws "should equal" getInterpolationWeightsBatch(grid, [1,2,3])
        @test_throws "size(x" getInterpolationWeightsBatch(grid, [1 1; 2 2; 3 3])
    end
    @testset "loadNeededPoints" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "dimension" loadNeededPoints!(grid, zeros(1))
        grid = makeSequenceGrid(dimension = 2, outputs = 2, depth = 2, type = "level", rule = "rleja")
        @test_throws "dimension" loadNeededPoints!(grid, zeros(3, 6))
        @test_throws "dimension" loadNeededPoints!(grid, zeros(2, 5))
        loadNeededPoints!(grid, zeros(2, 6))
        @test_throws "dimension" loadNeededPoints!(grid, ones(2,5))
        grid = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 1, order = 1, rule = "localp")
    end
    @testset "evaluate" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 2, depth = 2, type = "level", rule = "rleja")
        @test_throws "evaluate" evaluate(grid, zeros(2,1))
        @test_throws "evaluate" evaluateThreadSafe(grid, zeros(1, 2))
        @test_throws "evaluate" evaluateBatch(grid, zeros(2, 1))
        loadNeededPoints!(grid, zeros(2, 6))
        @test_throws "x should" evaluate(grid, zeros(3,1))
        @test_throws "x should" evaluate(grid, zeros(3))
        @test_throws "x should" evaluateThreadSafe(grid, zeros(1, 3))
        @test_throws "x should" evaluateThreadSafe(grid, zeros(3))
        @test_throws "vals" evaluateBatch(grid, zeros(3, 1))
        @test_throws "size(vals" evaluateBatch(grid, zeros(1, 2))
    end
    @testset "interpolate" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 2, depth = 2, type = "level", rule = "rleja")
        @test_throws "integrate" integrate(grid)
    end
    @testset "setDomainTransform!" begin
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "gauss-legendre")
        @test_throws "transformation" setDomainTransform!(grid, zeros(2))
        @test_throws "transformation" setDomainTransform!(grid, zeros(1,2))
        grid = makeGlobalGrid(dimension = 3, outputs = 1, depth = 2, type = "level", rule = "gauss-legendre")
        @test_throws "transformation" setDomainTransform!(grid, zeros(2,2))
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "gauss-legendre")
        @test_throws "transformation" setDomainTransform!(grid, zeros(3,2))
        grid = makeGlobalGrid(dimension = 3, outputs = 1, depth = 2, type = "level", rule = "clenshaw-curtis")
        @test setConformalTransformASIN!(grid, [0, 2, 4]) isa Nothing
        @test_throws "truncation" setConformalTransformASIN!(grid, [0, 2])
        @test_throws "truncation" setConformalTransformASIN!(grid, [0 2 3; 1 2 3])
    end
    @testset "setAnisotropicRefinement" begin
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "loadNeededPoints" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = 1)
        grid = makeGlobalGrid(dimension = 2, outputs = 0, depth = 2, type = "level", rule = "clenshaw-curtis");
        @test_throws "outputs" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = 1)
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2,  type = "level", rule = "fejer2");
        loadExpN2!(grid)
        @test_throws "min_growth" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = -2, output = 1)
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "clenshaw-curtis")
        loadExpN2!(grid)
        @test_throws "output" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = -1)
        loadExpN2!(grid)
        @test setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = 0) isa Nothing
        @test setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = 0, level_limits = [2, 3]) isa Nothing
        @test_throws "limits" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = 0, level_limits = [2, 3, 3])
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 3, type = "iptotal", rule = "leja")
        loadExpN2!(grid)
        @test setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = -1) isa Nothing
        @test_throws "output" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = -2)
        @test_throws "output" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = 5)
        @test_throws "type" setAnisotropicRefinement!(grid, type = "wrong", min_growth = 10, output = 0)
        @test setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = -1, level_limits = [3, 4]) isa Nothing
        @test_throws "limits" setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 10, output = -1,level_limits =  [3, 4, 5])
    end
    @testset "estimateAnisotropicCoefficients" begin
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
        @test_throws "loadNeededPoints" estimateAnisotropicCoefficients(grid, type = "iptotal", output = 1)
        grid = makeGlobalGrid(dimension = 2, outputs = 0, depth = 2, type = "level", rule = "clenshaw-curtis");
        @test_throws "outputs" estimateAnisotropicCoefficients(grid, type = "iptotal", output = 1)
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "clenshaw-curtis");
        loadExpN2!(grid);
        @test_throws "output" estimateAnisotropicCoefficients(grid, type = "iptotal", output = -1)
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 3, type = "iptotal", rule = "leja")
        loadExpN2!(grid);
        @test estimateAnisotropicCoefficients(grid, type = "iptotal", output = -1) isa Vector{Int32}
        @test estimateAnisotropicCoefficients(grid, type = "ipcurved", output = -1) isa Vector{Int32}
        @test_throws "output" estimateAnisotropicCoefficients(grid, type = "iptotal", output = -2)
        @test_throws "output" estimateAnisotropicCoefficients(grid, type = "iptotal", output = 5)
        @test_throws "type" estimateAnisotropicCoefficients(grid, type = "wrong", output = 0)
    end
    @testset "setSurplusRefinement!" begin
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "clenshaw-curtis")
        @test_throws "non-sequence" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic")
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "leja")
        @test_throws "loadNeededPoints" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic")
        loadExpN2!(grid)
        @test_throws "tolerance" setSurplusRefinement!(grid, tolerance = -1.E-4, output=0, refinement_type="classic")
        @test_throws "Sequence Grids" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic")
        @test  setSurplusRefinement!(grid, tolerance = 1.E-4, output=0) isa Nothing
        @test_throws "limits" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="", level_limits=[2, 3, 4])
        @test setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="", level_limits=[2, 3]) isa Nothing
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 2, order = 1, rule = "localp")
        loadExpN2!(grid)
        @test_throws "refinement_type" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0)
        @test_throws "refinement_type" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="class")
        @test setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic") isa Nothing
        @test_throws "limits" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic", level_limits=[2, 3, 4])
        @test_throws "scale" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic", level_limits=[], scale_correction=ones(1, 3))
        @test_throws "scale" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic", level_limits=[], scale_correction=ones(1, getNumPoints(grid) - 1))
        @test_throws "scale" setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic", level_limits=[], scale_correction=ones(2, getNumPoints(grid)))
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 2, depth = 2, order = 1, rule = "localp")
        loadExpN2!(grid)
        @test setSurplusRefinement!(grid, tolerance = 1.E-4, output=-1, refinement_type="classic", level_limits=[], scale_correction=ones(2, getNumPoints(grid))) isa Nothing
        @test_throws "scale" setSurplusRefinement!(grid, tolerance = 1.E-4, output=-1, refinement_type="classic", level_limits=[], scale_correction=ones(3, getNumPoints(grid)))
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 2, order = 1, rule = "localp")
        loadExpN2!(grid)
        @test setSurplusRefinement!(grid, tolerance = 1.E-4, output=0, refinement_type="classic", level_limits=[2, 3]) isa Nothing
    end
    @testset "removePointsByHierarchicalCoefficient" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "leja")
        @test_throws "polynomial" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = 0)
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 2, order = 1, rule = "localp")
        @test_throws "tolerance" removePointsByHierarchicalCoefficient!(grid, tolerance = -1.E-4, ouput = 0)
        @test_throws "output" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = -2)
        @test_throws "output" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = 3)
        @test_throws "loaded" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = 0)
        loadExpN2!(grid)
        @test removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = 0) isa Nothing
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 3, depth = 2, order = 1, rule = "localp")
        loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid)))
        @test_throws "scale" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = -1, scale_correction = ones(3))
        @test_throws "scale" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = -1,  scale_correction =ones(10))
        @test_throws "scale" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = -1,  scale_correction =ones(2,11))
        @test_throws "scale" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = -1,  scale_correction =ones(2,13))
        @test removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = -1,  scale_correction =ones(3,13)) isa Nothing
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 3, depth = 2, order = 1, rule = "localp")
        loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid)))
        @test_throws "scale" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = 0, scale_correction = ones(3,13))
        @test_throws "scale" removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = 0, scale_correction = ones(11))
        @test removePointsByHierarchicalCoefficient!(grid, tolerance = 1.E-4, output = 0, scale_correction = ones(1,13)) isa Nothing
    end
    @testset "evaluateHierarchicalFunctions" begin
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "clenshaw-curtis")
        @test evaluateHierarchicalFunctions(grid, [1.0 1.0; 0.5 0.3]) isa AbstractMatrix
        @test evaluateHierarchicalFunctions(grid, [1.0; 0.5]) isa Matrix
        @test evaluateHierarchicalFunctions(grid, [1.0, 1.0]) isa Matrix
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 2, order = 1, rule = "localp")
        @test evaluateSparseHierarchicalFunctions(grid, [1.0 1.0; 0.5 0.3]) isa AbstractMatrix
        @test_throws "getNumDimensions" evaluateSparseHierarchicalFunctions(grid, [1.0 0.5])
        @test_throws "matrix"  evaluateSparseHierarchicalFunctions(grid, [1.0, 1.0])
    end
    @testset "setHierarchicalCoefficients" begin
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 2, order = 1, rule = "localp")
        @test_throws "matrix" setHierarchicalCoefficients!(grid, [1.0, 1.0])
        @test_throws "columns" setHierarchicalCoefficients!(grid, [1.0 1.0])
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 1, order = 1, rule = "localp")
        @test_throws "columns" setHierarchicalCoefficients!(grid, [1.0 1.0; 1.0 1.0; 1.0 1.0; 1.0 1.0; 1.0 1.0])
        @test setHierarchicalCoefficients!(grid, ones(1, 5)) isa Nothing
        grid = makeFourierGrid(dimension = 2, outputs = 1, depth = 1, type = "level")
        @test_throws "Complex" setHierarchicalCoefficients!(grid, ones(1, 5))
    end
    @testset "Construction" begin
        grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 3, type = "level", rule = "rleja")
        @test_throws "before" getCandidateConstructionPoints(grid, type = "level", anisotropic_weights_or_output = -1)
        beginConstruction(grid)
        @test_throws "type" getCandidateConstructionPoints(grid, type = "lev", anisotropic_weights_or_output = -1)
        @test_throws "anisotropic_weights_or_output" getCandidateConstructionPoints(grid, type = "level", anisotropic_weights_or_output = [0])
        @test_throws "anisotropic_weights_or_output" getCandidateConstructionPoints(grid, type = "level", anisotropic_weights_or_output = "string")
        @test_throws "level_limits" getCandidateConstructionPoints(grid, type = "level", anisotropic_weights_or_output = -1, level_limits = [2])
        @test getCandidateConstructionPoints(grid, type = "level", anisotropic_weights_or_output = -1, level_limits = [2, 1]) isa Array
        grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 1, order = 1, rule = "localp")
        @test_throws "before" getCandidateConstructionPointsSurplus(grid, tolerance = 1.E-5, refinement_type = "classic")
        beginConstruction(grid)
        @test_throws "level_limits" getCandidateConstructionPointsSurplus(grid, tolerance = 1.E-5, refinement_type = "classic", output = 0, level_limit = [2])
        @test getCandidateConstructionPointsSurplus(grid, tolerance = 1.E-5, refinement_type = "classic", output = -1, level_limits = [2, 3]) isa Array{Float64}
        beginConstruction(grid);
        @test loadConstructedPoint!(grid, [0.0, 0.0], [1.0]) isa Nothing
        @test_throws "x" loadConstructedPoint!(grid, [0.0], [1.0])
        @test_throws "y" loadConstructedPoint!(grid, [0.0, 0.0], [1.0, 2.0])
        @test_throws "y" loadConstructedPoint!(grid, hcat([0.0, 0.0], [1.0, 0.0]), [1.0, 2.0])
        @test_throws "x" loadConstructedPoint!(grid, hcat([0.0,], [1.0,]), hcat([1.0,], [2.0,]))
        @test_throws "y" loadConstructedPoint!(grid, hcat([0.0, 0.0], [1.0, 0.0]), hcat([1.0, 2.0], [1.0, 2.0]))
        @test_throws "y" loadConstructedPoint!(grid, hcat([0.0, 0.0], [1.0, 0.0]), hcat([1.0,], [2.0,], [3.0,]))
        @test_throws "x" loadConstructedPoint!(grid, [hcat([0.0, 0.0], [1.0, 0.0])], hcat([1.0,], [2.0,]))
    end
    @testset "GPU" begin
        grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "leja")
        @test_throws "invalid acceleration" enableAcceleration!(grid, "gpu-wrong")
        @test enableAcceleration!(grid, "gpu-default") isa Nothing
        @test enableAcceleration!(grid, "gpu-default", GPUID=isAccelerationAvailable(grid, "gpu-cuda") ? 0 : nothing) isa Nothing
        @test_throws "GPU" enableAcceleration!(grid, "gpu-default", GPUID=-11)
        grid1 = TasmanianSG()
        @test_throws "invalid acceleration" isAccelerationAvailable(grid1, "cpu-wrong")
        @test isAccelerationAvailable(grid1, "cpu-blas") isa Bool
        @test_throws "GPU" getGPUMemory(grid1, -1)
        @test_throws "GPU" getGPUMemory(grid1, 1000000)
        @test_throws "GPU" getGPUMemory(grid1, getNumGPUs())
        @test_throws "GPU" getGPUName(grid1, -1)
        @test_throws "GPU" getGPUName(grid1, 1000000)
        @test_throws "GPU" getGPUName(grid1, getNumGPUs())
        @test_throws "GPU" setGPUID!(grid1, -1)
        @test_throws "GPU" setGPUID!(grid1, 1000000)
        @test_throws "GPU" setGPUID!(grid1, getNumGPUs())
    end
end

@testset verbose = true "Testing error handling" begin
   @testset "getSparseGridTests"  getSparseGridTests()
end
