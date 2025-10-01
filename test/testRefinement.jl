#=
import unittest
import TasmanianSG
import numpy as np

from random import shuffle

import testCommon

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
=#

"""
  Test the refinement capabilities:
        * set different refinements
        * estimate anisotropic coefficients
        * read/write refinements
"""
#=
def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")
=#
using Tasmanian
using Test

include("testCommon.jl")

"""
    Set refinement and clear refinement
"""
function checkSetClear()
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 2, order = 1, rule = "semi-localp")
    loadExpN2!(grid)
    @test getNumNeeded(grid) == 0
    setSurplusRefinement!(grid,  tolerance = 0.0001, output = 0, refinement_type = "classic")
    @test getNumNeeded(grid)  > 0 
    Tasmanian.clearRefinement!(grid) 
    @test getNumNeeded(grid) == 0
end

"""
        Check anisotropic coefficients
"""
function checkAnisoCoeff()
    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 9, type = "level", rule = "rleja")
    aP = getPoints(grid)
    aV = exp.(aP[1, :] + aP[2, :].^2)
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
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 4, order = 1, rule = "semi-localp")
    loadExpN2!(grid)
    aPoints = getPoints(grid)
    aScale = Matrix{Float64}((aPoints[1,:] .> 0)')
    setSurplusRefinement!(grid,  tolerance = 1.E-9, output = 0, refinement_type = "classic", level_limits = Int32[], scale_correction = aScale)
    aNeeded = getNeededPoints(grid)
    @test all(aNeeded[1, :] .> 0)
end

"""
        Read/Write regular refinement.
"""
function checkFileIO()
    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 2, type="level", rule="clenshaw-curtis")
    checkFileIO(grid)
    
    grid = makeSequenceGrid(dimension = 2, outputs = 1, depth = 2, type="level", rule="rleja")
    checkFileIO(grid)
end

function checkFileIO(grid)
    loadExpN2!(grid)
    setAnisotropicRefinement!(grid,  type = "iptotal", min_growth = 20, output = 0)
    @test getNumNeeded(grid)  > 0
    gridB = copyGrid(grid)
    compareGrids(grid, gridB)
    
    write(grid, "refTestFlename.grid", binary = true)
    gridB = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 0, order = 1)
    read(gridB, "refTestFlename.grid")
    compareGrids(grid, gridB)
    
    write(grid, "refTestFlename.grid", binary = false)
    gridB = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 0, order = 1)
    read(gridB, "refTestFlename.grid")
    compareGrids(grid, gridB)
 end   
    
#=
"""
        Test read/write when using construction.
"""
function checkConstruction()
        llTest = ["gridA.makeGlobalGrid(3, 2, 2, "level", "clenshaw-curtis"); gridB.makeGlobalGrid(3, 2, 2, "level", "clenshaw-curtis")",
                  "gridA.makeSequenceGrid(3, 2, 4, "level", "leja"); gridB.makeSequenceGrid(3, 2, 4, "level", "leja")",
                  "gridA.makeLocalPolynomialGrid(3, 2, 2); gridB.makeLocalPolynomialGrid(3, 2, 2)",
                  "gridA.makeWaveletGrid(3, 2, 2); gridB.makeWaveletGrid(3, 2, 2)",
                  "gridA.makeFourierGrid(3, 2, 2, "level"); gridB.makeFourierGrid(3, 2, 2, "level")",]

        for sMakeGrids in llTest:
            for sFormat in [False, True]: # test binary and ascii format
                gridA = TasmanianSG.TasmanianSparseGrid()
                gridB = TasmanianSG.TasmanianSparseGrid()
                gridC = TasmanianSG.TasmanianSparseGrid()

                exec(sMakeGrids)

                gridA.beginConstruction()
                gridB.beginConstruction()
                #gridA.printStats()

                gridB.write("testSave", bUseBinaryFormat = sFormat)
                gridB.makeSequenceGrid(1, 1, 0, "level", "rleja") # clean the grid
                gridB.read("testSave")
                compareGrids(gridA, gridB)
                gridC.copyGrid(gridA)
                compareGrids(gridA, gridC)

                for t in range(5): # use 5 iterations
                    if (gridA.isLocalPolynomial() or gridA.isWavelet()):
                        aPointsA = gridA.getCandidateConstructionPointsSurplus(1.E-4, "fds")
                        aPointsB = gridB.getCandidateConstructionPointsSurplus(1.E-4, "fds")
                        aPointsC = gridC.getCandidateConstructionPointsSurplus(1.E-4, "fds")
                    else:
                        aPointsA = gridA.getCandidateConstructionPoints("level", 0)
                        aPointsB = gridB.getCandidateConstructionPoints("level", 0)
                        aPointsC = gridC.getCandidateConstructionPoints("level", 0)
                    np.testing.assert_almost_equal(aPointsA, aPointsB, decimal=11)
                    np.testing.assert_almost_equal(aPointsA, aPointsC, decimal=11)

                    iNumPoints = int(aPointsA.shape[0] / 2)
                    if (iNumPoints > 32): iNumPoints = 32

                    # use the first samples (up to 32) and shuffle the order
                    # add one of the samples further in the list
                    liSamples = list(range(iNumPoints + 1))
                    shuffle(liSamples)
                    for iI in range(len(liSamples)):
                        if (liSamples[iI] == iNumPoints):
                            liSamples[iI] = iNumPoints + 1
                    #liSamples = map(lambda i: i if i < iNumPoints else iNumPoints + 1, liSamples)

                    for iI in liSamples: # compute and load the samples
                        aPoint = aPointsA[iI, :]
                        aValue = np.array([exp(aPoint[0] + aPoint[1]), 1.0 / ((aPoint[0] - 1.3) * (aPoint[1] - 1.6) * (aPoint[2] - 2.0))])

                        gridA.loadConstructedPoint(aPoint, aValue)
                        gridB.loadConstructedPoint(aPoint, aValue)
                        gridC.loadConstructedPoint(aPoint, aValue)

                    # using straight construction or read/write should produce the same result
                    compareGrids(gridA, gridC)
                    gridB.write("testSave", bUseBinaryFormat = sFormat)
                    gridB.makeSequenceGrid(1, 1, 0, "level", "rleja")
                    gridB.read("testSave")
                    compareGrids(gridA, gridB)
                    gridC.copyGrid(gridA)
                    compareGrids(gridA, gridC)

                gridA.finishConstruction()
                gridB.finishConstruction()

                gridB.write("testSave", bUseBinaryFormat = sFormat)
                gridB.makeSequenceGrid(1, 1, 0, "level", "rleja")
                gridB.read("testSave")
                compareGrids(gridA, gridB)
                gridC.copyGrid(gridA)
                compareGrids(gridA, gridC)

        # check multi-point load
        gridA = TasmanianSG.TasmanianSparseGrid()
        gridA.makeLocalPolynomialGrid(3, 2, 4);
        loadExpN2!(gridA)

        gridB = TasmanianSG.TasmanianSparseGrid()
        gridB.makeLocalPolynomialGrid(3, 2, 0)

        gridB.beginConstruction()
        aX = gridA.getPoints()
        aY = gridA.evaluateBatch(aX)
        gridB.loadConstructedPoint(aX, aY)
        gridB.finishConstruction()
        compareGrids(gridA, gridB)

        # check some mem-leaks and crashes (correctness is elsewhere)
        gridA = TasmanianSG.TasmanianSparseGrid()
        gridA.makeLocalPolynomialGrid(2, 5, 0)
        gridA.beginConstruction()
        gridA.loadConstructedPoint(np.empty([0, 2]), np.empty([0, 5])) # empty input, check for crash

        gridA.makeLocalPolynomialGrid(2, 1, 1)
        gridA.loadNeededPoints(ones([5, 1]))
        gridA.beginConstruction()
        aPoints = gridA.getCandidateConstructionPointsSurplus(1.E-4, "classic") # should generate empty output
        np.testing.assert_almost_equal(aPoints, np.empty([0, 0]), 14, "failed to generate empty list of construction points", True)

        gridA.makeLocalPolynomialGrid(2, 1, 0)
        gridA.loadNeededPoints(ones([1, 1]))
        gridA.beginConstruction()
        aPoints = gridA.getCandidateConstructionPointsSurplus(1.E-4, "classic", 0, [], np.array([[1.E-6]])) # should generate empty output
        np.testing.assert_almost_equal(aPoints, np.empty([0, 0]), 14, "failed to generate empty list of construction points", True)

        gridA.makeGlobalGrid(2, 1, 1, "tensor", "clenshaw-curtis")
        gridA.loadNeededPoints(ones([9, 1]))
        gridA.beginConstruction()
        aPoints = gridA.getCandidateConstructionPoints("ipcurved", [5, 5, 2, 2], [1, 1]) # should generate empty output
        np.testing.assert_almost_equal(aPoints, np.empty([0, 0]), 14, "failed to generate empty list of construction points", True)
=#

"""
        tests removePointsByHierarchicalCoefficient()
"""
function checkRemovePoints()
    grid = makeLocalPolynomialGrid(dimension = 2, outputs = 1, depth = 1)
    aPoints = getNeededPoints(grid) 
    loadNeededValues!(grid, exp.(-aPoints[1, :].^2 -0.5*aPoints[2, :].^2))

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
    removePointsByHierarchicalCoefficient!(reduced, tolerance = 0.0, output = -1, scale_correction = [], num_new_points = 3)
    @test getNumPoints(reduced) == 3
    @test getLoadedPoints(reduced) ≈ hcat([0.0, 0.0], [-1.0, 0.0], [1.0, 0.0])

    reduced = copyGrid(grid)
    removePointsByHierarchicalCoefficient!(reduced, tolerance = 0.0, output = -1, scale_correction = [], num_new_points = 1)
    @test getNumPoints(reduced) ==  1
    @test getLoadedPoints(reduced) ≈ [0.0, 0.0]

    reduced = copyGrid(grid)
    removePointsByHierarchicalCoefficient!(reduced, tolerance = 0.0, output = -1, scale_correction = hcat([1.0], [1.0], [1.0], [0.1], [0.1]), num_new_points = 3)
    @test getNumPoints(reduced) == 3
    @test getLoadedPoints(reduced) ≈ hcat([0.0, 0.0], [0.0, -1.0], [0.0, 1.0])
end


@testset verbose=true "Testing core refine grid" begin
    @testset "Set refinement and clear refinement" checkSetClear()
    
    @testset "Check anisotropic coefficients" checkAnisoCoeff()
    
    @testset "Check surplus refinement for local polynomial grids" checkLocalpSurplus()
    
    @testset "Read/Write regular refinement" checkFileIO()
    
    #=
    @testset "Test Read/Write when using construction" checkConstruction()
    =#

    @testset "Tests removePointsByHierarchicalCoefficient()" checkRemovePoints()
end

nothing    
