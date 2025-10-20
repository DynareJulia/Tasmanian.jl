import Random: shuffle!
using Tasmanian
using Test

"""
    Set refinement and clear refinement
"""
function checkSetClear()
    grid = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 2, order = 1, rule = "semi-localp")
    loadExpN2!(grid)
    @test getNumNeeded(grid) == 0
    setSurplusRefinement!(grid, tolerance = 0.0001, output = 0, refinement_type = "classic")
    @test getNumNeeded(grid) > 0
    Tasmanian.clearRefinement!(grid)
    @test getNumNeeded(grid) == 0
end

"""
        Check anisotropic coefficients
"""
function checkAnisoCoeff()
    grid = makeGlobalGrid(
        dimension = 2, outputs = 1, depth = 9, type = "level", rule = "rleja")
    aP = getPoints(grid)
    aV = exp.(aP[1, :] + aP[2, :] .^ 2)
    loadNeededPoints!(grid, aV)
    aC = estimateAnisotropicCoefficients(grid, type = "iptotal", output = 0)
    @test ndims(aC) == 1
    @test size(aC, 1) == 2
    @test abs(aC[1] / aC[2] - 2.0) < 0.2
    aC = estimateAnisotropicCoefficients(grid, type = "ipcurved", output = 0)
    @test ndims(aC) == 1
    @test size(aC, 1) == 4
    @test abs(aC[1] / aC[2] - 2.0) < 0.2
    @test aC[3] < 0.0
    @test aC[4] < 0.0
end

"""
        Check surplus refinement for local polynomial grids
"""
function checkLocalpSurplus()
    grid = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 4, order = 1, rule = "semi-localp")
    loadExpN2!(grid)
    aPoints = getPoints(grid)
    aScale = Matrix{Float64}((aPoints[1, :] .> 0)')
    setSurplusRefinement!(grid, tolerance = 1.E-9, output = 0, refinement_type = "classic",
        level_limits = Int32[], scale_correction = aScale)
    aNeeded = getNeededPoints(grid)
    @test all(aNeeded[1, :] .> 0)
end

"""
        Read/Write regular refinement.
"""
function checkFileIO()
    grid = makeGlobalGrid(
        dimension = 2, outputs = 1, depth = 2, type = "level", rule = "clenshaw-curtis")
    checkFileIO(grid)

    grid = makeSequenceGrid(
        dimension = 2, outputs = 1, depth = 2, type = "level", rule = "rleja")
    checkFileIO(grid)
end

function checkFileIO(grid)
    loadExpN2!(grid)
    setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 20, output = 0)
    @test getNumNeeded(grid) > 0
    gridB = copyGrid(grid)
    compareGrids(grid, gridB)

    write(grid, "refTestFlename.grid", binary = true)
    gridB = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 0, order = 1)
    read!(gridB, "refTestFlename.grid")
    compareGrids(grid, gridB)

    write(grid, "refTestFlename.grid", binary = false)
    gridB = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 0, order = 1)
    read!(gridB, "refTestFlename.grid")
    compareGrids(grid, gridB)
end

function checkConstruction_(gridA, gridB)
    for format in [false, true] # test binary and ascii format
        gridC = TasmanianSG()

        beginConstruction!(gridA)
        beginConstruction!(gridB)

        write(gridB, "testSave", binary = format)
        makeSequenceGrid!(
            gridB, dimension = 1, outputs = 1, depth = 0, type = "level", rule = "rleja") # clean the grid
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)
        copyGrid!(gridC, gridA)
        compareGrids(gridA, gridC)

        for t in 1:5 # use 5 iterations
            if isLocalPolynomial(gridA) || isWavelet(gridA)
                PointsA = getCandidateConstructionPointsSurplus(
                    gridA, tolerance = 1.E-4, refinement_type = "fds")
                PointsB = getCandidateConstructionPointsSurplus(
                    gridB, tolerance = 1.E-4, refinement_type = "fds")
                PointsC = getCandidateConstructionPointsSurplus(
                    gridC, tolerance = 1.E-4, refinement_type = "fds")
            else
                PointsA = getCandidateConstructionPoints(
                    gridA, type = "level", anisotropic_weights_or_output = 0)
                PointsB = getCandidateConstructionPoints(
                    gridB, type = "level", anisotropic_weights_or_output = 0)
                PointsC = getCandidateConstructionPoints(
                    gridC, type = "level", anisotropic_weights_or_output = 0)
            end
            @test PointsA ≈ PointsB
            @test PointsA ≈ PointsC

            NumPoints = Int(floor(size(PointsA, 2) / 2))
            NumPoints > 32 && (NumPoints = 32)

            # use the first samples (up to 32) and shuffle the order
            # add one of the samples further in the list
            samples = collect(1:NumPoints)
            shuffle!(samples)
            for iI in 1:length(samples)
                samples[iI] == NumPoints && (samples[iI] = NumPoints + 1)
            end

            for iI in samples # compute and load the samples
                Point = PointsA[:, iI]
                Value = [exp(Point[1] + Point[2]),
                    1.0 / ((Point[1] - 1.3) * (Point[2] - 1.6) * (Point[3] - 2.0))]

                loadConstructedPoint!(gridA, x = Point, y = Value)
                loadConstructedPoint!(gridB, x = Point, y = Value)
                loadConstructedPoint!(gridC, x = Point, y = Value)
            end

            # using straight construction or read/write should produce the same result
            compareGrids(gridA, gridC)
            write(gridB, "testSave", binary = format)
            makeSequenceGrid!(gridB, dimension = 1, outputs = 1,
                depth = 0, type = "level", rule = "rleja")
            read!(gridB, "testSave")
            compareGrids(gridA, gridB)
            copyGrid!(gridC, gridA)
            compareGrids(gridA, gridC)
        end

        finishConstruction!(gridA)
        finishConstruction!(gridB)

        write(gridB, "testSave", binary = format)
        makeSequenceGrid!(
            gridB, dimension = 1, outputs = 1, depth = 0, type = "level", rule = "rleja")
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)
        copyGrid!(gridC, gridA)
        compareGrids(gridA, gridC)
    end
end

"""
        Test read/write when using construction.
"""
function checkConstruction()
    gridA = TasmanianSG()
    gridB = TasmanianSG()

    makeGlobalGrid!(gridA, dimension = 3, outputs = 2, depth = 2,
        type = "level", rule = "clenshaw-curtis")
    makeGlobalGrid!(gridB, dimension = 3, outputs = 2, depth = 2,
        type = "level", rule = "clenshaw-curtis")
    checkConstruction_(gridA, gridB)

    makeSequenceGrid!(
        gridA, dimension = 3, outputs = 2, depth = 4, type = "level", rule = "leja")
    makeSequenceGrid!(
        gridB, dimension = 3, outputs = 2, depth = 4, type = "level", rule = "leja")
    checkConstruction_(gridA, gridB)

    makeLocalPolynomialGrid!(gridA, dimension = 3, outputs = 2, depth = 2)
    makeLocalPolynomialGrid!(gridB, dimension = 3, outputs = 2, depth = 2)
    checkConstruction_(gridA, gridB)

    makeWaveletGrid!(gridA, dimension = 3, outputs = 2, depth = 2)
    makeWaveletGrid!(gridB, dimension = 3, outputs = 2, depth = 2)
    checkConstruction_(gridA, gridB)

    makeFourierGrid!(gridA, dimension = 3, outputs = 2, depth = 2, type = "level")
    makeFourierGrid!(gridB, dimension = 3, outputs = 2, depth = 2, type = "level")
    checkConstruction_(gridA, gridB)

    # check multi-point load
    makeLocalPolynomialGrid!(gridA, dimension = 3, outputs = 2, depth = 4)
    loadExpN2!(gridA)

    makeLocalPolynomialGrid!(gridB, dimension = 3, outputs = 2, depth = 0)

    beginConstruction!(gridB)
    X = getPoints(gridA)
    Y = evaluateBatch(gridA, X)
    loadConstructedPoint!(gridB, x = X, y = Y)
    finishConstruction!(gridB)
    compareGrids(gridA, gridB)

    # check some mem-leaks and crashes (correctness is elsewhere)
    makeLocalPolynomialGrid!(gridA, dimension = 2, outputs = 5, depth = 0)
    beginConstruction!(gridA)
    loadConstructedPoint!(gridA, x = zeros(2, 0), y = zeros(5, 0)) # empty input, check for crash

    makeLocalPolynomialGrid!(gridA, dimension = 2, outputs = 1, depth = 1)
    loadNeededPoints!(gridA, ones(1, 5))
    beginConstruction!(gridA)
    Points = getCandidateConstructionPointsSurplus(
        gridA, tolerance = 1.E-4, refinement_type = "classic") # should generate empty output
    @test Points ≈ zeros(2, 0)

    makeLocalPolynomialGrid!(gridA, dimension = 2, outputs = 1, depth = 0)
    loadNeededPoints!(gridA, ones(1, 1))
    beginConstruction!(gridA)
    Points = getCandidateConstructionPointsSurplus(
        gridA, tolerance = 1.E-4, refinement_type = "classic",
        output = 0, scale_correction = [1.E-6]) # should generate empty output
    @test Points ≈ zeros(2, 0)

    makeGlobalGrid!(gridA, dimension = 2, outputs = 1, depth = 1,
        type = "tensor", rule = "clenshaw-curtis")
    loadNeededPoints!(gridA, ones(1, 9))
    beginConstruction!(gridA)
    Points = getCandidateConstructionPoints(gridA, type = "ipcurved",
        anisotropic_weights_or_output = [5, 5, 2, 2], level_limits = [1, 1]) # should generate empty output
    @test Points ≈ zeros(2, 0)
end

"""
        tests removePointsByHierarchicalCoefficient()
"""
function checkRemovePoints()
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 1)
    aPoints = getNeededPoints(grid)
    loadNeededValues!(grid, exp.(-aPoints[1, :] .^ 2 - 0.5 * aPoints[2, :] .^ 2))

    reduced = TasmanianSG(2, 1, 1)
    reduced = copyGrid(grid)
    removePointsByHierarchicalCoefficient!(reduced, tolerance = 0.6)
    @test getNumPoints(reduced) == 3
    @test getLoadedPoints(reduced) ≈ hcat([0.0, 0.0], [-1.0, 0.0], [1.0, 0.0])

    reduced = copyGrid(grid)
    removePointsByHierarchicalCoefficient!(reduced, tolerance = 0.7)
    @test getNumPoints(reduced) == 1
    @test getLoadedPoints(reduced) ≈ [0.0, 0.0]

    reduced = copyGrid(grid)
    removePointsByHierarchicalCoefficient!(
        reduced, tolerance = 0.0, output = -1, scale_correction = [], num_new_points = 3)
    @test getNumPoints(reduced) == 3
    @test getLoadedPoints(reduced) ≈ hcat([0.0, 0.0], [-1.0, 0.0], [1.0, 0.0])

    reduced = copyGrid(grid)
    removePointsByHierarchicalCoefficient!(
        reduced, tolerance = 0.0, output = -1, scale_correction = [], num_new_points = 1)
    @test getNumPoints(reduced) == 1
    @test getLoadedPoints(reduced) ≈ [0.0, 0.0]

    reduced = copyGrid(grid)
    removePointsByHierarchicalCoefficient!(reduced, tolerance = 0.0, output = -1,
        scale_correction = hcat([1.0], [1.0], [1.0], [0.1], [0.1]), num_new_points = 3)
    @test getNumPoints(reduced) == 3
    @test getLoadedPoints(reduced) ≈ hcat([0.0, 0.0], [0.0, -1.0], [0.0, 1.0])
end

@testset verbose=true "Testing core refine grid" begin
    @testset "Set refinement and clear refinement" checkSetClear()

    @testset "Check anisotropic coefficients" checkAnisoCoeff()

    @testset "Check surplus refinement for local polynomial grids" checkLocalpSurplus()

    @testset "Read/Write regular refinement" checkFileIO()

    @testset "Test Read/Write when using construction" checkConstruction()

    @testset "Tests removePointsByHierarchicalCoefficient()" checkRemovePoints()
end

nothing
