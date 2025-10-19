using Random
using Tasmanian, Tasmanian_jll
using Test

include("testCommon.jl")
include("testConfigureData.jl")

#=
    Test the acceleration options:
    * consistency between CMake options and options reported by C++/Python
    * ability to switch between GPUs, read properties set GPUID, etc.
    * consistency between accelerated and reference evaluations
=#

"""
    Check if the CMake options are propagated to the C++ library and the
    correct acceleration options are passed between C++ and Julia.
"""
function checkMeta()
    grid = TasmanianSG()

    if bEnableSyncTests
        @test isAccelerationAvailable(grid, "cpu-blas") == bHasBlas
        @test isAccelerationAvailable(grid, "gpu-cublas") == bHasCuBlas
        @test isAccelerationAvailable(grid, "gpu-cuda") == bHasCuda
        @test isAccelerationAvailable(grid, "gpu-default") == (bHasCuBlas || bHasCuda)

        # acceleration meta-data
        lsAvailableAcc = []
        bHasBlas && push!(lsAvailableAcc, "cpu-blas")
        bHasCuBlas && push!(lsAvailableAcc, "gpu-cublas")
        bHasCuda && push!(lsAvailableAcc, "gpu-cuda")
        grid = makeLocalPolynomialGrid(
            dimension = 2, outputs = 1, depth = 2, order = 1, rule = "semi-localp")
        for accel in lsAvailableAcc
            enableAcceleration!(grid, accel)
            sA = getAccelerationType(grid)
            @test accel == sA
        end
        lsAvailableAcc = []
        bHasCuBlas && push!(lsAvailableAcc, ("gpu-rocblas", "gpu-cublas"))
        bHasCuda && push!(lsAvailableAcc, ("gpu-hip", "gpu-cuda"))
        for accel in lsAvailableAcc
            enableAcceleration!(grid, accel[1])
            sA = getAccelerationType(grid)
            @test accel[2] == sA
        end
    end
end

"""
    Check setting and resetting the GPU ID and reading the device names.
"""
function checkMultiGPU()
    grid = TasmanianSG()

    if (bHasSycl)
        @test getGPUID(grid) == -1
    else
        @test getGPUID(grid) == 0

        if getNumGPUs() > 1
            setGPUID(grid, 1)
            @test getGPUID(grid) == 1
            if !bUsingMSVC
                sName = getGPUName(grid, 1) # mostly checks for memory leaks and crashes
            end
        end
        if getNumGPUs() > 0
            setGPUID(grid, 0)
            @test getGPUID(grid) == 0
            if !bUsingMSVC
                sName = getGPUName(grid, 0) # mostly checks for memory leaks and crashes
            end
        end
    end
end

"""
Check for consistency between accelerated and reference evaluations.
In short, set a grid, load points, evaluate in different ways.
Test all visible GPUs and combinations cuBlas/MAGMA, etc.
"""
function checkEvaluateConsistency()
    grid = TasmanianSG()

    aTestPointsCanonical = 2 * rand(2, 100) .- 1
    aTestPointsTransformed = 2 * rand(2, 100) .+ 3
    aDomainTransform = [3.0 5.0; 3.0 5.0]

    grid = makeGlobalGrid(
        dimension = 2, outputs = 2, depth = 4, type = "level", rule = "clenshaw-curtis")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeGlobalGrid(
        dimension = 2, outputs = 2, depth = 4, type = "level", rule = "chebyshev")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeSequenceGrid(
        dimension = 2, outputs = 2, depth = 4, type = "level", rule = "leja")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeLocalPolynomialGrid(
        dimension = 2, outputs = 3, depth = 4, order = 1, rule = "localp")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeLocalPolynomialGrid(
        dimension = 2, outputs = 2, depth = 4, order = 2, rule = "localp")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 4, order = 3, rule = "localp")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 4, order = 4, rule = "semi-localp")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeWaveletGrid(dimension = 2, outputs = 1, depth = 3, order = 1)
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeWaveletGrid(dimension = 2, outputs = 1, depth = 3, order = 3)
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    grid = makeFourierGrid(dimension = 2, outputs = 1, depth = 3, type = "level")
    checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
end

function checkEvaluateConsistency(
        grid, aTestPointsCanonical, aTestPointsTransformed, aDomainTransform)
    iNumGPUs = getNumGPUs()
    lsAccelTypes = ["none", "cpu-blas", "gpu-cuda", "gpu-cublas", "gpu-magma"]
    iFastEvalSubtest = 6

    for iI in 1:2
        iC = 1
        iGPU = 0
        iGPUID > -1 && (iGPU = iGPUID)

        while iC < length(lsAccelTypes)
            sAcc = lsAccelTypes[iC]
            if iI == 1
                aTestPoints = aTestPointsCanonical
            else
                aTestPoints = aTestPointsTransformed
                setDomainTransform!(grid, aDomainTransform)
            end
            enableAcceleration!(grid, sAcc)
            if (sAcc == "gpu-cuda") || (sAcc == "gpu-cublas") || (sAcc == "gpu-magma")
                if iGPU < getNumGPUs() # without cuda or cublas, NumGPUs is 0 and cannot set GPU
                    setGPUID!(grid, iGPU)
                end
            end

            loadExpN2!(grid)

            aRegular = stack(evaluateThreadSafe(grid, aTestPoints[:, i])
            for i in axes(aTestPoints, 2))
            aBatched = evaluateBatch(grid, aTestPoints)
            @test aRegular ≈ aBatched

            aFast = stack(evaluate(grid, aTestPoints[:, i]) for i in 1:iFastEvalSubtest)
            @test aRegular[:, 1:iFastEvalSubtest] ≈ aFast

            if (sAcc == "gpu-cuda") || (sAcc == "gpu-cublas") || (sAcc == "gpu-magma")
                if iGPUID == -1
                    iGPU += 1
                    if iGPU >= iNumGPUs
                        iC += 1
                        iGPU = 0
                    end
                else
                    iC += 1
                end
            else
                iC += 1
            end
        end
    end
end

@testset verbose=true "Testing accelerated evaluate consistency" begin
    @testset "checkMeta" checkMeta()

    @testset "checkMultiGPU" checkMultiGPU()

    @testset "checkEvaluateConsistency" checkEvaluateConsistency()
end
