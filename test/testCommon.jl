using Test

"""
        compareGrids(gridA, gridB; bTestRuleNames = true)

Compares two grids by checking points, weights, and several evaluate.
The evaluates are done on a canonical [-1, 1] interval.

The test passes if the grids are mathematically identical.

bTestRuleNames should be true in most cases, i.e., check if the names
of the rule in each grid matches. The exception is when comparing the
GaussPatterson grids build with custom file and sRule = "gauss-patterson"
"""
function compareGrids(gridA, gridB; bTestRuleNames = true)
    # test basic number of points data
    @assert getNumDimensions(gridA) == getNumDimensions(gridB) "error in getNumDimensions()"
    @assert getNumOutputs(gridA) == getNumOutputs(gridB) "error in getNumOutputs()"
    @assert getNumPoints(gridA) == getNumPoints(gridB) "error in getNumPoints()"
    @assert getNumLoaded(gridA) == getNumLoaded(gridB) "error in getNumLoaded()"
    @assert getNumNeeded(gridA) == getNumNeeded(gridB) "error in getNumNeeded()"
    @assert isUsingConstruction(gridA) == isUsingConstruction(gridB) "error in isUsingConstruction()"
    
    if getNumPoints(gridA) == 0 # emptry grid, nothing else to check
        return
    end
    # load test points (canonical domain only), make sure to avoid grid points, e.g., 1/2, 1/4, etc.
    if getNumDimensions(gridA) == 1
        mX1 = [1.0/3.0]
        mX2 = [-1.0/3.0]
        mX3 = [-1.0/5.0]
        mX4 = [1.0/7.0]
        mX5 = [-1.0/7.0]
    elseif getNumDimensions(gridA) == 2
        mX1 = [1.0/3.0, 1.0/6.0]
        mX2 = [-1.0/3.0, 1.0/6.0]
        mX3 = [-1.0/5.0, -1.0/7.0]
        mX4 = [1.0/7.0, 1.0/5.0]
        mX5 = [-1.0/7.0, -1.0/13.0]
    elseif getNumDimensions(gridA) == 3
        mX1 = [1.0/3.0, 1.0/6.0, -1.0/5.0]
        mX2 = [-1.0/3.0, 1.0/6.0, 1.0/6.0]
        mX3 = [-1.0/5.0, -1.0/7.0, 1.0/3.0]
        mX4 = [1.0/7.0, 1.0/5.0, 2.0/3.0]
        mX5 = [-1.0/7.0, -1.0/13.0, -2.0/3.0]
    end
    aBatchPoints = hcat(mX1, mX2, mX3, mX4, mX5)
    
    pA = getPoints(gridA)
    pB = getPoints(gridB)
    @assert pA == pB "Points not equal"
    
    pA = getLoadedPoints(gridA)
    pB = getLoadedPoints(gridB)
    @assert pA == pB "Loaded not equal"
    
    pA = getNeededPoints(gridA)
    pB = getNeededPoints(gridB)
    @assert pA == pB "Needed not equal"
    
    # test rule data
    @assert getAlpha(gridA) == getAlpha(gridB) "error in getAlpha()"
    @assert getBeta(gridA) == getBeta(gridB) "error in getBeta()"
    @assert getOrder(gridA) == getOrder(gridB) "error in getOrder()"

    if (bTestRuleNames)
        @assert getRule(gridA) == getRule(gridB) "error in getRule()"
        @assert getCustomRuleDescription(gridA) == getCustomRuleDescription(gridB) "error in getCustomRuleDescription()"
    end
    @assert isGlobal(gridA) == isGlobal(gridB) "error in isGlobal()"
    @assert isSequence(gridA) == isSequence(gridB) "error in isSequence()"
    @assert isLocalPolynomial(gridA) == isLocalPolynomial(gridB) "error in isLocalPolynomial()"
    @assert isWavelet(gridA) == isWavelet(gridB) "error in isWavelet()"
    @assert isFourier(gridA) == isFourier(gridB) "error in isFourier()"

    # weights
    pA = getQuadratureWeights(gridA)
    pB = getQuadratureWeights(gridB)
    @assert pA == pB "Quadrature not equal"

    pA = getInterpolationWeights(gridA, mX1)
    pB = getInterpolationWeights(gridB, mX1)
    @assert pA == pB "Interpolation test 1 not equal"
    
    pA = getInterpolationWeights(gridA, mX2)
    pB = getInterpolationWeights(gridB, mX2)
    @assert pA== pB "Interpolation test 2 not equal"
    
    pA = getInterpolationWeights(gridA, mX3)
    pB = getInterpolationWeights(gridB, mX3)
    @assert pA == pB "Interpolation test 3 not equal"

    # evaluate (values have been loaded)
    if (getNumLoaded(gridA) > 0)
        pA = evaluate(gridA, mX4)
        pB = evaluate(gridB, mX4)
        @assert pA ≈ pB "Interpolation test 4 not equal"
        
        pA = evaluate(gridA, mX5)
        pB = evaluate(gridB, mX5)
        @assert pA ≈ pB "Interpolation test 5 not equal"

        pA = integrate(gridA)
        pB = integrate(gridB)
        @assert pA ≈ pB "Integration test not equal"

        pA = evaluateBatch(gridA, aBatchPoints)
        pB = evaluateBatch(gridB, aBatchPoints)
        @assert pA ≈ pB "Interpolation test 6 (batch) not equal"

        pA = getHierarchicalCoefficients(gridA)
        pB = getHierarchicalCoefficients(gridB)
        @assert pA ≈ pB "getHierarchicalCoefficients() not equal"

    end

    # domain transforms
    @assert isSetDomainTransform(gridA) == isSetDomainTransform(gridB) "error in isSetDomainTransfrom()"
    
    pA = getDomainTransform(gridA)
    pB = getDomainTransform(gridB)
    @assert pA == pB "Domain test no equal"
    
    @assert isSetConformalTransformASIN(gridA) == isSetConformalTransformASIN(gridB) "error in isSetConformalTransformASIN()"
    
    pA = getLevelLimits(gridA)
    pB = getLevelLimits(gridB)
    @assert pA == pB "Level limit test no equal"
    
    pA = getConformalTransformASIN(gridA)
    pB = getConformalTransformASIN(gridB)
    @assert pA == pB "Conformal transform ASIN not equal"

    return true
end

import Base.==
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
    @assert getDescription(ctA) == getDescription(ctB) "error in getDescription()"
    @assert getNumLevels(ctA) == getNumLevels(ctB) "error in getNumLevels()"
    for level in 1:ctA.getNumLevels()
        @assert ctA.getNumPoints(level) == ctB.getNumPoints(level) "error in getNumPoints() at level $level"
        @assert ctA.getIExact(level) == ctB.getIExact(level) "error in getIExact() at level $level"
        @assert ctA.getQExact(level) == ctB.getQExact(level) "error in getQExact() at level $level"
        wA, nA = ctA.getWeightsNodes(level)
        wB, nB = ctB.getWeightsNodes(level)
        @assert wA == wB "Weights from getWeightsNodes() are not equal at level $level" 
        @assert nA == nB "Nodes from getWeightsNodes() are not equal at level $level"
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
    loadNeededPoints!(grid, repeat(exp.(-sum(mPoints.^2, dims=1)), iOut))
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
