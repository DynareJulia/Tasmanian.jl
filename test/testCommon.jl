import Base.==
using Test

"""
    compareGrids(gridA, gridB; bTestRuleNames = true)

Compares two grids by checking points, weights, and several evaluate.
The evaluates are done on a canonical [-1, 1] interval.

The test passes if the grids are mathematically identical.

bTestRuleNames should be true in most cases, i.e., check if the names
of the rule in each grid matches. The exception is when comparing the
GaussPatterson grids build with custom file and rule = "gauss-patterson"
"""
function compareGrids(gridA::TasmanianSG, gridB::TasmanianSG; bTestRuleNames = true)
    # test basic number of points data
    @test getNumDimensions(gridA) == getNumDimensions(gridB)
    @test getNumOutputs(gridA) == getNumOutputs(gridB)
    @test getNumPoints(gridA) == getNumPoints(gridB)
    @test getNumLoaded(gridA) == getNumLoaded(gridB)
    @test getNumNeeded(gridA) == getNumNeeded(gridB)
    @test isUsingConstruction(gridA) == isUsingConstruction(gridB)

    if getNumPoints(gridA) == 0 # emptry grid, nothing else to check
        return
    end
    # load test points (canonical domain only), make sure to avoid grid points, e.g., 1/2, 1/4, etc.
    if getNumDimensions(gridA) == 1
        mX1 = [1.0 / 3.0]
        mX2 = [-1.0 / 3.0]
        mX3 = [-1.0 / 5.0]
        mX4 = [1.0 / 7.0]
        mX5 = [-1.0 / 7.0]
    elseif getNumDimensions(gridA) == 2
        mX1 = [1.0 / 3.0, 1.0 / 6.0]
        mX2 = [-1.0 / 3.0, 1.0 / 6.0]
        mX3 = [-1.0 / 5.0, -1.0 / 7.0]
        mX4 = [1.0 / 7.0, 1.0 / 5.0]
        mX5 = [-1.0 / 7.0, -1.0 / 13.0]
    elseif getNumDimensions(gridA) == 3
        mX1 = [1.0 / 3.0, 1.0 / 6.0, -1.0 / 5.0]
        mX2 = [-1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0]
        mX3 = [-1.0 / 5.0, -1.0 / 7.0, 1.0 / 3.0]
        mX4 = [1.0 / 7.0, 1.0 / 5.0, 2.0 / 3.0]
        mX5 = [-1.0 / 7.0, -1.0 / 13.0, -2.0 / 3.0]
    else
        mX1 = mX2 = mX3 = mX4 = mX5 = Float64[]
    end
    aBatchPoints = hcat(mX1, mX2, mX3, mX4, mX5)

    pA = getPoints(gridA)
    pB = getPoints(gridB)
    @test pA == pB

    pA = getLoadedPoints(gridA)
    pB = getLoadedPoints(gridB)
    @test pA == pB

    pA = getNeededPoints(gridA)
    pB = getNeededPoints(gridB)
    @test pA == pB

    # test rule data
    @test getAlpha(gridA) == getAlpha(gridB)
    @test getBeta(gridA) == getBeta(gridB)
    @test getOrder(gridA) == getOrder(gridB)

    if (bTestRuleNames)
        @test getRule(gridA) == getRule(gridB)
        @test getCustomRuleDescription(gridA) == getCustomRuleDescription(gridB)
    end
    @test isGlobal(gridA) == isGlobal(gridB)
    @test isSequence(gridA) == isSequence(gridB)
    @test isLocalPolynomial(gridA) == isLocalPolynomial(gridB)
    @test isWavelet(gridA) == isWavelet(gridB)
    @test isFourier(gridA) == isFourier(gridB)

    # weights
    pA = getQuadratureWeights(gridA)
    pB = getQuadratureWeights(gridB)
    @test pA == pB

    pA = getInterpolationWeights(gridA, mX1)
    pB = getInterpolationWeights(gridB, mX1)
    @test pA == pB

    pA = getInterpolationWeights(gridA, mX2)
    pB = getInterpolationWeights(gridB, mX2)
    @test pA == pB

    pA = getInterpolationWeights(gridA, mX3)
    pB = getInterpolationWeights(gridB, mX3)
    @test pA == pB

    # evaluate (values have been loaded)
    if (getNumLoaded(gridA) > 0)
        pA = evaluate(gridA, mX4)
        pB = evaluate(gridB, mX4)
        @test pA ≈ pB

        pA = evaluate(gridA, mX5)
        pB = evaluate(gridB, mX5)
        @test pA ≈ pB

        pA = integrate(gridA)
        pB = integrate(gridB)
        @test pA ≈ pB

        pA = evaluateBatch(gridA, aBatchPoints)
        pB = evaluateBatch(gridB, aBatchPoints)
        @test pA ≈ pB

        pA = getHierarchicalCoefficients(gridA)
        pB = getHierarchicalCoefficients(gridB)
        @test pA ≈ pB
    end

    # domain transforms
    @test isSetDomainTransform(gridA) == isSetDomainTransform(gridB)

    pA = getDomainTransform(gridA)
    pB = getDomainTransform(gridB)
    @test pA == pB

    @test isSetConformalTransformASIN(gridA) == isSetConformalTransformASIN(gridB)

    pA = getLevelLimits(gridA)
    pB = getLevelLimits(gridB)
    @test pA == pB

    pA = getConformalTransformASIN(gridA)
    pB = getConformalTransformASIN(gridB)
    @test pA == pB

    return true
end

"""
 ==(a::TasmanianSG, b::TasmaninanSG)

is shorthand for compareGrids(gridA, gridB)
"""
==(gridA::TasmanianSG, gridB::TasmanianSG) = compareGrids(gridA, gridB)

"""
    Compares two CustomTabulated instances by checking nodes, weights, and other metadata.
    
The test passes if the two instances are mathematically identical.
"""
function compareCustomTabulated(ctA, ctB)
    # Test metadata.
    @test getDescription(ctA) == getDescription(ctB)
    @test getNumLevels(ctA) == getNumLevels(ctB)
    for level in 0:(getNumLevels(ctA) - 1)
        @test getNumPoints(ctA, level) == getNumPoints(ctB, level)
        @test getIExact(ctA, level) == getIExact(ctB, level)
        @test getQExact(ctA, level) == getQExact(ctB, level)
        wA, nA = getWeightsNodes(ctA, level)
        wB, nB = getWeightsNodes(ctB, level)
        @test wA == wB
        @test nA == nB
    end
end

"""
If there are needed points, load the points with exp(-sum(x)**2) where
x is each of the needed points.

The function is needed to test whether read/write correctly works on the
loaded values, so some values have to be easily loaded even if they
are not meaningful as in the convergence/correctness tests.
"""
function loadExpN2!(grid)
    if getNumNeeded(grid) == 0
        return
    end
    mPoints = getNeededPoints(grid)
    iOut = getNumOutputs(grid)
    loadNeededPoints!(grid, repeat(exp.(-sum(mPoints .^ 2, dims = 1)), iOut))
end

"""
aPoints is 1D array of points
MustHave and MustNotHave are the points to test
Example: checkPoints(aPoints[:,0], [0.0, 1.0, -1.0], [0.5, -0.5])
"""
function checkPoints(aPoints; MustHave, MustNotHave)
    for x in MustHave
        @test any(abs.(aPoints .- x) .< 0.001)
    end

    for x in MustNotHave
        @test !any(abs.(aPoints .- x) .< 0.001)
    end
end
