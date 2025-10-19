using Tasmanian, Tasmanian_jll
using Test

include("testCommon.jl")

"""
    Test make/update, check points, check library paths
"""

"""
    Check different ways to specify the library path and check the
    version number.
"""
function checkPathsVersions()
    grid = TasmanianSG()
    v = pkgversion(Tasmanian_jll)
    @test Tasmanian.getVersionMajor() == v.major
    @test Tasmanian.getVersionMinor() == v.minor
    sLicense = Tasmanian.getLicense()
    LicenseFile = read(
        joinpath(Tasmanian_jll.artifact_dir, "share", "licenses", "Tasmanian", "LICENSE"),
        String)
    @test sLicense[1:12] == LicenseFile[1:12]
end

"""
    Make grids and check against pen-and-paper answers, check weights too.
"""
function checkMakeAgainstKnown()
    grid = makeGlobalGrid(dimension = 1, outputs = 0, depth = 4, type = "level",
        rule = "gauss-hermite", anisotropic_weights = [], alpha = 2.0)
    aW = getQuadratureWeights(grid)
    @test abs(sum(aW) - 0.5 * sqrt(pi)) < 1.E-14

    grid = makeGlobalGrid(dimension = 2, outputs = 0, depth = 2, type = "level",
        rule = "leja", anisotropic_weights = [2, 1])
    aA = hcat([0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0])
    aP = getPoints(grid)
    @test aA == aP

    grid = makeGlobalGrid(dimension = 2, outputs = 0, depth = 4, type = "ipcurved",
        rule = "leja", anisotropic_weights = [20, 10, 0, 7])
    aA = hcat([0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [0.0, sqrt(1.0 / 3.0)],
        [1.0, 0.0], [1.0, 1.0], [-1.0, 0.0])
    aP = getPoints(grid)
    @test aA ≈ aP

    grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level",
        rule = "leja", anisotropic_weights = [2, 1])
    aA = hcat([0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0])
    aP = getPoints(grid)
    @test aA == aP

    grid = Tasmanian.makeSequenceGrid(
        dimension = 2, outputs = 1, depth = 4, type = "ipcurved",
        rule = "leja", anisotropic_weights = [20, 10, 0, 7])
    aA = hcat([0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [0.0, sqrt(1.0 / 3.0)],
        [1.0, 0.0], [1.0, 1.0], [-1.0, 0.0])
    aP = getPoints(grid)
    @test aA ≈ aP

    grid = Tasmanian.makeFourierGrid(
        dimension = 2, outputs = 1, depth = 2, type = "level", anisotropic_weights = [1, 2])
    aA = hcat([0.0, 0.0], [0.0, 1.0 / 3.0], [0.0, 2.0 / 3.0], [1.0 / 3.0, 0.0],
        [2.0 / 3.0, 0.0], [1.0 / 9.0, 0.0], [2.0 / 9.0, 0.0],
        [4.0 / 9.0, 0.0], [5.0 / 9.0, 0.0], [7.0 / 9.0, 0.0], [8.0 / 9.0, 0.0])
    aP = getPoints(grid)
    @test aA ≈ aP

    # this is a very important test, checks the curved rule and covers the non-lower-set index selection code
    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 1, type = "ipcurved",
        rule = "rleja", anisotropic_weights = [10, 10, -21, -21])
    aA = hcat([0.0, 0.0], [0.0, 1.0], [0.0, 0.5], [0.0, 0.25], [0.0, 0.75], [0.0, 0.125],
        [1.0, 0.0], [1.0, 1.0], [1.0, 0.5], [1.0, 0.25], [1.0, 0.75], [1.0, 0.125],
        [0.5, 0.0], [0.5, 1.0], [0.5, 0.5], [0.5, 0.25], [0.5, 0.75], [0.5, 0.125],
        [0.25, 0.0], [0.25, 1.0], [0.25, 0.5], [0.25, 0.25], [0.25, 0.75],
        [0.75, 0.0], [0.75, 1.0], [0.75, 0.5], [0.75, 0.25],
        [0.125, 0.0], [0.125, 1.0], [0.125, 0.5])
    aA = cos.(pi * aA)
    aP = getPoints(grid)
    @test aA ≈ aP

    # test weights
    aX = hcat([0.0, -0.3], [-0.44, 0.7], [0.82, -0.01])
    grid = makeGlobalGrid(
        dimension = 2, outputs = 0, depth = 4, type = "level", rule = "chebyshev")
    aW = getInterpolationWeightsBatch(grid, aX)
    aX1 = aX[:, 1]
    aX2 = aX[:, 2]
    aX3 = aX[:, 3]
    aW1 = getInterpolationWeights(grid, aX1)
    aW2 = getInterpolationWeights(grid, aX2)
    aW3 = getInterpolationWeights(grid, aX3)
    aWW = hcat(aW1, aW2, aW3)
    @test aW == aWW

    grid = makeGlobalGrid(
        dimension = 3, outputs = 0, depth = 1, type = "level", rule = "chebyshev")
    aA = hcat([0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
    @test aA == getPoints(grid)
    aTrans = vcat([-1.0 -4.0], [2.0 5.0], [-3.0 5])
    setDomainTransform!(grid, aTrans)
    aB = hcat([-2.5, 3.5, 1.0], [-2.5, 3.5, -3.0], [-2.5, 3.5, 5.0],
        [-2.5, 2.0, 1.0], [-2.5, 5.0, 1.0], [-1.0, 3.5, 1.0], [-4.0, 3.5, 1.0])
    @test aB == getPoints(grid)
    clearDomainTransform!(grid)
    @test aA == getPoints(grid)

    grid = makeGlobalGrid(
        dimension = 3, outputs = 0, depth = 1, type = "level", rule = "fejer2")
    #display(getPoints(grid))
    aA = hcat(
        [0.0, 0.0, 0.0], [0.0, 0.0, -0.707106781186548], [0.0, 0.0, 0.707106781186548],
        [0.0, -0.707106781186548, 0.0], [0.0, 0.707106781186548, 0.0],
        [-0.707106781186548, 0.0, 0.0], [0.707106781186548, 0.0, 0.0])
    #display(aA)
    @test aA ≈ getPoints(grid)
    setConformalTransformASIN!(grid, Int32[4, 6, 3])
    aA = hcat([0.0, 0.0, 0.0], [0.0, 0.0, -0.60890205], [0.0, 0.0, 0.60890205],
        [0.0, -0.57892643, 0.0], [0.0, 0.57892643, 0.0],
        [-0.59587172, 0.0, 0.0], [0.59587172, 0.0, 0.0])
    @test aA ≈ getPoints(grid)
    clearConformalTransform!(grid)
    aA = hcat(
        [0.0, 0.0, 0.0], [0.0, 0.0, -0.707106781186548], [0.0, 0.0, 0.707106781186548],
        [0.0, -0.707106781186548, 0.0], [0.0, 0.707106781186548, 0.0],
        [-0.707106781186548, 0.0, 0.0], [0.707106781186548, 0.0, 0.0])
    @test aA ≈ getPoints(grid)

    # number of points
    grid = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 2, order = 2, rule = "localp")
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 13
    @test iNL == 0
    @test iNP == iNN
    loadExpN2!(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 0
    @test iNL == 13
    @test iNP == iNL
    setSurplusRefinement!(grid, tolerance = 0.10, output = 0, refinement_type = "classic")
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 8
    @test iNL == 13
    @test iNP == iNL
    setSurplusRefinement!(grid, tolerance = 0.05, output = 0, refinement_type = "classic")
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 16
    @test iNL == 13
    @test iNP == iNL
    loadExpN2!(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 0
    @test iNL == 29
    @test iNP == iNL
end

"""
    Check the update and make working together.
"""
function checkUpdate()
    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level",
        rule = "leja", anisotropic_weights = [2, 1])
    checkUpdate_1(grid)
    updateGlobalGrid!(grid, depth = 2, type = "level", anisotropic_weights = [2, 1])
    checkUpdate_2(grid)
    updateGlobalGrid!(grid, depth = 3, type = "level", anisotropic_weights = [2, 1])
    checkUpdate_3(grid)
    updateGlobalGrid!(
        grid, depth = 4, type = "ipcurved", anisotropic_weights = [20, 10, 0, 7])
    checkUpdate_4(grid)

    grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type = "level",
        rule = "leja", anisotropic_weights = [2, 1])
    checkUpdate_1(grid)
    updateSequenceGrid!(grid, depth = 2, type = "level", anisotropic_weights = [2, 1])
    checkUpdate_2(grid)
    updateSequenceGrid!(grid, depth = 3, type = "level", anisotropic_weights = [2, 1])
    checkUpdate_3(grid)
    updateSequenceGrid!(
        grid, depth = 4, type = "ipcurved", anisotropic_weights = [20, 10, 0, 7])
    checkUpdate_4(grid)

    grid = makeFourierGrid(dimension = 2, outputs = 1, depth = 2, type = "level")
    loadExpN2!(grid)
    updateFourierGrid!(grid, depth = 3, type = "level")
    @test getNumNeeded(grid) == 60
end

function checkUpdate_1(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 4
    @test iNL == 0
    @test iNP == iNN
    loadExpN2!(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 0
    @test iNL == 4
    @test iNP == iNL
end

function checkUpdate_2(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 0
    @test iNL == 4
    @test iNP == iNL
end

function checkUpdate_3(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 2
    @test iNL == 4
    @test iNP == iNL
    aA = hcat([0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0])
    aP = getPoints(grid)
    @test aA == aP
    aP = getLoadedPoints(grid)
    @test aA == aP
    aA = hcat([0.0, sqrt(1.0 / 3.0)], [1.0, 1.0])
    aP = getNeededPoints(grid)
    @test aA == aP
    loadExpN2!(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 0
    @test iNL == 6
    @test iNP == iNL
end

function checkUpdate_4(grid)
    iNN = getNumNeeded(grid)
    iNL = getNumLoaded(grid)
    iNP = getNumPoints(grid)
    @test iNN == 1
    @test iNL == 6
    @test iNP == iNL
    @test iNN == 1
    @test iNL == 6
    @test iNP == iNL
    aA = hcat([-1.0, 0.0])
    aP = getNeededPoints(grid)
    @test aA == aP
end

"""
    Test that default values are set correctly for alpha/beta/order,
    empty returns, empty copy.
"""
function checkDefaults()
    grid = TasmanianSG()
    for sRule in Tasmanian.LocalRules
        grid = makeLocalPolynomialGrid(
            dimension = 3, outputs = 1, depth = 3, order = 0, rule = sRule)
        @test getAlpha(grid) == 0.0
        @test getBeta(grid) == 0.0
        @test getOrder(grid) == 0
        grid = makeLocalPolynomialGrid(
            dimension = 3, outputs = 1, depth = 3, order = 1, rule = sRule)
        @test getAlpha(grid) == 0.0
        @test getBeta(grid) == 0.0
        @test getOrder(grid) == 1
        grid = makeLocalPolynomialGrid(
            dimension = 3, outputs = 1, depth = 3, order = 2, rule = sRule)
        @test getAlpha(grid) == 0.0
        @test getBeta(grid) == 0.0
        @test getOrder(grid) == 2
    end
    grid = makeWaveletGrid(dimension = 3, outputs = 1, depth = 2, order = 1)
    @test getAlpha(grid) == 0.0
    @test getBeta(grid) == 0.0
    @test getOrder(grid) == 1
    grid = makeWaveletGrid(dimension = 3, outputs = 1, depth = 2, order = 3)
    @test getAlpha(grid) == 0.0
    @test getBeta(grid) == 0.0
    @test getOrder(grid) == 3

    # default alpha/beta and order
    grid = makeGlobalGrid(dimension = 2, outputs = 0, depth = 2, type = "level",
        rule = "leja", anisotropic_weights = [2, 1])
    @test getAlpha(grid) == 0.0
    @test getBeta(grid) == 0.0
    @test getOrder(grid) == -1
    grid = makeGlobalGrid(dimension = 1, outputs = 0, depth = 4,
        type = "level", rule = "gauss-hermite", alpha = 2.0)
    @test getAlpha(grid) == 2.0
    @test getBeta(grid) == 0.0
    @test getOrder(grid) == -1
    grid = makeGlobalGrid(dimension = 1, outputs = 0, depth = 4, type = "level",
        rule = "gauss-jacobi", alpha = 3.0, beta = 2.0)
    @test getAlpha(grid) == 3.0
    @test getBeta(grid) == 2.0
    @test getOrder(grid) == -1

    grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 1, type = "level",
        rule = "leja", anisotropic_weights = [2, 1])
    @test getAlpha(grid) == 0.0
    @test getBeta(grid) == 0.0
    @test getOrder(grid) == -1

    # test empty returns
    dummy_ans = zeros(0, 0)
    grid_dummy = TasmanianSG()
    aX = zeros(0)
    @test dummy_ans == getPoints(grid_dummy)
    @test dummy_ans == getNeededPoints(grid_dummy)
    @test dummy_ans == getLoadedPoints(grid_dummy)
    dummy_ans = zeros(0)
    @test dummy_ans == getQuadratureWeights(grid_dummy)
    #@test dummy_ans == getInterpolationWeights(grid_dummy, aX)
    dummy_ans = zeros(0, 0)
    aX = zeros(0, 0)
    #@test dummy_ans == getInterpolationWeightsBatch(grid_dummy, aX)
    aX = zeros(1, 0)
    #@test dummy_ans == getInterpolationWeightsBatch(grid_dummy, aX)
    grid = makeGlobalGrid(
        dimension = 2, outputs = 1, depth = 2, type = "level", rule = "chebyshev")
    loadExpN2!(grid)
    dummy_ans = zeros(1, 0)
    aX = zeros(2, 0)
    @test dummy_ans == evaluateBatch(grid, aX)
end

"""
    Check setting level limits.
"""
function checkLevelLimits()
    grid = TasmanianSG()

    # Level Limits
    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type = "level",
        rule = "clenshaw-curtis", level_limits = [1, 4])
    aPoints = getPoints(grid)
    checkPoints(aPoints[1, :], MustHave = [-1.0, 1.0],
        MustNotHave = [1.0 / sqrt(2.0), -1.0 / sqrt(2.0)])
    checkPoints(aPoints[2, :], MustHave = [-1.0, 1.0, 1.0 / sqrt(2.0), -1.0 / sqrt(2.0)],
        MustNotHave = [])
    loadExpN2!(grid)
    setAnisotropicRefinement!(
        grid, type = "iptotal", min_growth = 20, output = 0, level_limits = [1, 4])
    @test getNumNeeded(grid) > 0
    aPoints = getNeededPoints(grid)
    checkPoints(
        aPoints[1, :], MustHave = [], MustNotHave = [1.0 / sqrt(2.0), -1.0 / sqrt(2.0)])

    grid = makeSequenceGrid(dimension = 3, outputs = 1, depth = 3, type = "level",
        rule = "leja", level_limits = [3, 2, 1])
    aPoints = getPoints(grid)
    checkPoints(
        aPoints[1, :], MustHave = [0.0, -1.0, 1.0, 1.0 / sqrt(3.0)], MustNotHave = [])
    checkPoints(aPoints[2, :], MustHave = [0.0, -1.0, 1.0], MustNotHave = [1.0 / sqrt(3.0)])
    checkPoints(aPoints[3, :], MustHave = [0.0, 1.0], MustNotHave = [-1.0, 1.0 / sqrt(3.0)])
    loadExpN2!(grid)
    setAnisotropicRefinement!(grid, type = "iptotal", min_growth = 5, output = 0)
    @test getNumNeeded(grid) > 0
    aPoints = getNeededPoints(grid)
    checkPoints(aPoints[2, :], MustHave = [], MustNotHave = [1.0 / sqrt(3.0)])
    checkPoints(aPoints[3, :], MustHave = [], MustNotHave = [-1.0, 1.0 / sqrt(3.0)])
    clearRefinement!(grid)
    setAnisotropicRefinement!(
        grid, type = "iptotal", min_growth = 10, output = 0, level_limits = [3, 2, 2])
    @test getNumNeeded(grid) > 0
    aPoints = getNeededPoints(grid)
    checkPoints(aPoints[2, :], MustHave = [], MustNotHave = [1.0 / sqrt(3.0)])
    checkPoints(aPoints[3, :], MustHave = [1.0], MustNotHave = [1.0 / sqrt(3.0)])
    clearRefinement!(grid)
    setSurplusRefinement!(
        grid, tolerance = 1.E-8, output = 0, refinement_type = "", level_limits = [3, 2, 1])
    @test getNumNeeded(grid) > 0
    aPoints = getNeededPoints(grid)
    checkPoints(aPoints[2, :], MustHave = [], MustNotHave = [1.0 / sqrt(3.0)])
    checkPoints(aPoints[3, :], MustHave = [], MustNotHave = [-1.0, 1.0 / sqrt(3.0)])

    # check that nodes from level 2 (+-0.5) and 3 (+-0.25, +-0.75) appear only in the proper dimension
    grid = makeLocalPolynomialGrid(dimension = 3, outputs = 1, depth = 3, order = 1,
        rule = "localp", level_limits = [1, 2, 3])
    aPoints = getPoints(grid)
    checkPoints(aPoints[1, :], MustHave = [0.0, -1.0, 1.0],
        MustNotHave = [0.5, -0.5, -0.75, -0.25, 0.25, 0.75])
    checkPoints(aPoints[2, :], MustHave = [0.0, -1.0, 1.0, 0.5, -0.5],
        MustNotHave = [-0.75, -0.25, 0.25, 0.75])
    checkPoints(
        aPoints[3, :], MustHave = [0.0, -1.0, 1.0, 0.5, -0.5, -0.75, -0.25, 0.25, 0.75],
        MustNotHave = [])
    loadExpN2!(grid)
    setSurplusRefinement!(grid, tolerance = 1.E-8, output = 0, refinement_type = "classic")
    @test getNumNeeded(grid) > 0
    aPoints = getNeededPoints(grid)
    checkPoints(
        aPoints[1, :], MustHave = [], MustNotHave = [0.5, -0.5, -0.75, -0.25, 0.25, 0.75])
    checkPoints(aPoints[2, :], MustHave = [], MustNotHave = [-0.75, -0.25, 0.25, 0.75])
    checkPoints(aPoints[3, :], MustHave = [], MustNotHave = [])
    clearRefinement!(grid)
    setSurplusRefinement!(grid, tolerance = 1.E-8, output = 0,
        refinement_type = "classic", level_limits = [2, 2, 3])
    aPoints = getNeededPoints(grid)
    checkPoints(
        aPoints[1, :], MustHave = [0.5, -0.5], MustNotHave = [-0.75, -0.25, 0.25, 0.75])
    checkPoints(aPoints[2, :], MustHave = [], MustNotHave = [-0.75, -0.25, 0.25, 0.75])

    grid = makeWaveletGrid(
        dimension = 2, outputs = 1, depth = 3, order = 1, level_limits = [0, 2])
    aPoints = getPoints(grid)
    checkPoints(aPoints[1, :], MustHave = [0.0, -1.0, 1.0], MustNotHave = [-0.5, 0.5])
    checkPoints(
        aPoints[2, :], MustHave = [0.0, -1.0, 1.0, -0.5, 0.5, -0.75, -0.25, 0.25, 0.75],
        MustNotHave = [-0.125, 0.125])
    loadExpN2!(grid)
    setSurplusRefinement!(grid, tolerance = 1.E-8, output = 0, refinement_type = "classic") # grid is effectively full tensor, cannot refine within the level limits
    @test getNumNeeded(grid) == 0
    setSurplusRefinement!(grid, tolerance = 1.E-8, output = 0,
        refinement_type = "classic", level_limits = [1, 2])
    aPoints = getNeededPoints(grid)
    checkPoints(
        aPoints[1, :], MustHave = [0.5, -0.5], MustNotHave = [-0.75, -0.25, 0.25, 0.75])
    checkPoints(aPoints[2, :], MustHave = [], MustNotHave = [-0.125, 0.125])

    # level limits I/O
    grid = makeLocalPolynomialGrid(dimension = 3, outputs = 1, depth = 3, order = 1,
        rule = "localp", level_limits = [1, 2, 3])
    aLimits = getLevelLimits(grid)
    @test aLimits ≈ [1, 2, 3]

    clearLevelLimits!(grid)
    aLimits = getLevelLimits(grid)
    @test aLimits ≈ [-1, -1, -1]
end

@testset verbose=true "Testing core make/update grid" begin
    @testset "checkPathsVersions" checkPathsVersions()

    @testset "checkMakeAgainstKnown" checkMakeAgainstKnown()

    @testset "checkUpdate" checkUpdate()

    @testset "checkDefaults" checkDefaults()

    @testset "checkLevelLimits" checkLevelLimits()
end
