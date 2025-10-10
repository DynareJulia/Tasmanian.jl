using LinearAlgebra
using Tasmanian, Tasmanian_jll
using Test

include("testCommon.jl")

        """
        Check hierarchical functions against known (pen-and-paper)
        solutions.
        """
function checkAgainsKnown()
    # evalHierarchicalBasis (all done in batch), dense mode
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 2, depth = 1, order = 1, rule = "localp")
    aT = zeros(2,0)
    @test aT == getHierarchicalCoefficients(grid)
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 0, depth = 1, order = 1, rule = "localp")
    aT = zeros(0,0)
    @test aT == getHierarchicalCoefficients(grid)
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 2, depth = 1, order = 1, rule = "localp")
    aPoints = getPoints(grid)
    aV = vcat(exp.(2.0 * aPoints[1,:] + aPoints[2,:])', sin.(3.0 * aPoints[1,:] + aPoints[2,:])')
    loadNeededPoints!(grid, aV)
    aS = getHierarchicalCoefficients(grid)
    aT = hcat([1.0, 0.0], [exp(-1.0)-1.0, sin(-1.0)], [exp(1.0)-1.0, sin(1.0)], [exp(-2.0)-1.0, sin(-3.0)], [exp(2.0)-1.0, sin(3.0)])
    @test aS ≈ aT

    grid = makeGlobalGrid(dimension = 3, outputs = 1, depth = 4, type = "level", rule = "fejer2")
    aPoints = getPoints(grid)
    aVan = evaluateHierarchicalFunctions(grid, aPoints)
    @test aVan ≈ I

    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 1, type = "level", rule = "clenshaw-curtis")
    aPoints = hcat([0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666])
    f0 = (1.0 .- aPoints[1,:].^2) + (1.0 .- aPoints[2,:].^2) .- 1.0
    f1 = 0.5 * aPoints[2,:] .* (aPoints[2,:] .- 1.0)
    f2 = 0.5 * aPoints[2,:] .* (aPoints[2,:] .+ 1.0)
    f3 = 0.5 * aPoints[1,:] .* (aPoints[1,:] .- 1.0)
    f4 = 0.5 * aPoints[1,:] .* (aPoints[1,:] .+ 1.0)
    aResult = vcat(f0', f1', f2', f3', f4')
    aVan = evaluateHierarchicalFunctions(grid, aPoints)
    @test aVan ≈ aResult

    grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level", rule = "leja")
    aPoints = hcat([0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666])
    f0 = ones(size(aPoints, 2))
    f1 = aPoints[2,:]
    f2 = 0.5 * aPoints[2,:] .* (aPoints[2,:] .- 1.0)
    f3 = aPoints[1,:]
    f4 = aPoints[1,:] .* aPoints[2,:]
    f5 = 0.5 * aPoints[1,:] .* (aPoints[1,:] .- 1.0)
    aResult = vcat(f0', f1', f2', f3', f4', f5')
    aVan = evaluateHierarchicalFunctions(grid, aPoints)
    @test aVan ≈ aResult

    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 1, order = 1, rule = "localp")
    aPoints = hcat([0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666])

    f0 = ones(size(aPoints, 2))
    f1 = -aPoints[2,:]
    f2 = copy(aPoints[2,:])
    f3 = -aPoints[1,:]
    f4 = copy(aPoints[1,:])
    for iI in axes(aPoints, 2)
        (aPoints[2,iI] > 0.0) && (f1[iI] = 0.0)
        (aPoints[2,iI] < 0.0) && (f2[iI] = 0.0)
        (aPoints[1,iI] > 0.0) && (f3[iI] = 0.0)
        (aPoints[1,iI] < 0.0) && (f4[iI] = 0.0)
    end
    aResult = vcat(f0', f1', f2', f3', f4')
    aVan = evaluateHierarchicalFunctions(grid, aPoints)
    @test aVan ≈ aResult

    grid = makeFourierGrid(dimension = 2, outputs = 1, depth = 1, type = "level")
    aPoints = [0.25, 0.5]
    aResult = [1.0, -1.0, -1.0, -1im, 1im]
    aVan = evaluateHierarchicalFunctions(grid, aPoints)
    @test aVan ≈ aResult

    # sparse hierarchical functions
    spmat = Tasmanian.getDenseForm(Tasmanian.TasmanianSimpleSparseMatrix())
    @test spmat == zeros(0,0)

    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 4, order = 1, rule = "localp")
    aPoints = hcat([0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666])
    aDense = evaluateHierarchicalFunctions(grid, aPoints)
    pSparse = evaluateSparseHierarchicalFunctions(grid, aPoints)
    @test aDense ≈ Matrix(pSparse)

    grid = makeWaveletGrid(dimension = 2, outputs = 1, depth = 4, order = 1)
    aPoints = hcat([0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666])
    aDense = evaluateHierarchicalFunctions(grid, aPoints)
    pSparse = evaluateSparseHierarchicalFunctions(grid, aPoints)
    @test aDense ≈ Matrix(pSparse)
end

"""
    Set the coefficients and use mergeRefinement()
"""
function checkSetCoeffsMergeRegine()
    aPoints = hcat([0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666])
    # set hierarchical coeffs/merge refinement
    gridA = makeGlobalGrid(dimension = 2, outputs = 1, depth = 4, type = "level", rule = "chebyshev")
    gridB = makeGlobalGrid(dimension = 2, outputs = 1, depth = 4, type = "level", rule = "chebyshev")
    loadExpN2!(gridA)
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)

    gridA = makeGlobalGrid(dimension = 2, outputs = 1, depth = 4, type = "level", rule = "fejer2")
    gridB = makeGlobalGrid(dimension = 2, outputs = 1, depth = 4, type = "level", rule = "fejer2")
    loadExpN2!(gridA)
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)
    setAnisotropicRefinement!(gridA, type = "iptotal", min_growth = 10, output = 0)
    setAnisotropicRefinement!(gridB, type = "iptotal", min_growth = 10, output = 0)
    loadExpN2!(gridA)
    mergeRefinement!(gridB)
    @test getNumPoints(gridA) == getNumPoints(gridB)
    aRes = evaluateBatch(gridB, aPoints)
    @test aRes ≈ zeros(size(aRes))
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)

    gridA = makeSequenceGrid(dimension = 2, outputs = 1, depth = 5, type = "ipcurved", rule = "min-delta")
    gridB = makeSequenceGrid(dimension = 2, outputs = 1, depth = 5, type = "ipcurved", rule = "min-delta")
    loadExpN2!(gridA)
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)
    setAnisotropicRefinement!(gridA, type = "iptotal", min_growth = 30, output = 0)
    setAnisotropicRefinement!(gridB, type = "iptotal", min_growth = 30, output = 0)
    loadExpN2!(gridA)
    mergeRefinement!(gridB)
    @test getNumPoints(gridA) == getNumPoints(gridB)
    aRes = evaluateBatch(gridB, aPoints)
    @test aRes ≈ zeros(size(aRes))
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)

    gridA = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 4, order = 2, rule = "semi-localp")
    gridB = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 4, order = 2, rule = "semi-localp")
    loadExpN2!(gridA)
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)
    setSurplusRefinement!(gridA, tolerance = 1.E-4, output = 0, refinement_type = "classic")
    setSurplusRefinement!(gridB, tolerance = 1.E-4, output = 0, refinement_type = "classic")
    loadExpN2!(gridA)
    mergeRefinement!(gridB)
    @test getNumPoints(gridA) == getNumPoints(gridB)
    aRes = evaluateBatch(gridB, aPoints)
    @test aRes ≈ zeros(size(aRes))
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)

    gridA = makeWaveletGrid(dimension = 2, outputs = 1, depth = 2, order = 1)
    gridB = makeWaveletGrid(dimension = 2, outputs = 1, depth = 2, order = 1)
    loadExpN2!(gridA)
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)
    setSurplusRefinement!(gridA, tolerance = 1.E-3, output = 0, refinement_type = "classic")
    setSurplusRefinement!(gridB, tolerance = 1.E-3, output = 0, refinement_type = "classic")
    loadExpN2!(gridA)
    mergeRefinement!(gridB)
    @test getNumPoints(gridA) == getNumPoints(gridB)
    aRes = evaluateBatch(gridB, aPoints)
    @test aRes ≈ zeros(size(aRes))
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)

    gridA = makeFourierGrid(dimension = 2, outputs = 1, depth = 3, type = "level")
    gridB = makeFourierGrid(dimension = 2, outputs = 1, depth = 3, type = "level")
    loadExpN2!(gridA)
    setHierarchicalCoefficients!(gridB, getHierarchicalCoefficients(gridA))
    compareGrids(gridA, gridB)
end

function checkIntegrals()
    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 4, type = "level", rule = "clenshaw-curtis")
    checkIntegrals(grid)
    grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 4, type = "level", rule = "rleja")
    checkIntegrals(grid)
    grid = makeFourierGrid(dimension = 2, outputs = 1, depth = 3, type = "level")
    checkIntegrals(grid)
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 4)
    checkIntegrals(grid)
    grid = makeWaveletGrid(dimension = 2, outputs = 1, depth = 2)
    checkIntegrals(grid)
end

function checkIntegrals(grid)
    setDomainTransform!(grid, vcat([-3.0 5.0], [-4.0 2.0]))
    loadExpN2!(grid)
    aIntegrals = integrateHierarchicalFunctions(grid)
    aSurps = getHierarchicalCoefficients(grid)
    fResult = aSurps * aIntegrals
    fExpected = integrate(grid)
    @test fResult[1] ≈ fExpected[1]

    # check support
    grid = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 3)
    aSup = getHierarchicalSupport(grid)
    aResult = [1.0 1.0 1.0 0.5 0.5 0.25 0.25 0.25 0.25]
    @test aSup ≈ aResult
end

function checkValues()
    grid = TasmanianSG()
    aResult = getLoadedValues(grid)
    @test length(aResult) == 0
    grid = makeGlobalGrid(dimension = 2, outputs = 2, depth = 3, type = "level", rule = "clenshaw-curtis")
    grid = makeSequenceGrid(dimension = 2, outputs = 2, depth = 3, type = "level", rule = "rleja")
    grid = makeFourierGrid(dimension = 2, outputs = 2, depth = 2, type = "level")
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 2, depth = 2)
    grid = makeWaveletGrid(dimension = 2, outputs = 2, depth = 2)
end

function checkValues_(grid)
    iNP = getNumNeeded(grid)
    iOuts = 2
    aReferece = [float(i) for i in 1:iNP * iOuts]
    loadNeededPoints(grid, aReferece)
    aResult = getLoadedValues(grid)
    @test aResult ≈ aReferece
end

@testset verbose = true "Testing core learning from random samples" begin
    @testset "checkAgainsKnown" checkAgainsKnown()

    @testset "checkSetCoeffsMergeRefine" checkSetCoeffsMergeRegine()

    @testset "checIntegrals" checkIntegrals()

    @testset "checkValues" checkValues()
end
