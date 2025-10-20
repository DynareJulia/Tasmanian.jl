using Format
using Tasmanian
using Test

"""
    Test the version I/O, available acceleration, gpu info, print stats
"""
function checkMetaIO()
    grid = TasmanianSG()
    format("Tasmanian Sparse Grids version: {1:1s}", Tasmanian.getVersion())
    format("                 version major: {1:1d}", Tasmanian.getVersionMajor())
    format("                 version minor: {1:1d}", Tasmanian.getVersionMinor())
    format("                       License: {1:1s}", Tasmanian.getLicense())
    if Tasmanian.isOpenMPEnabled()
        status = "Enabled"
    else
        status = "Disabled"
    end
    format("                        OpenMP: {1:1s}", status)
    AvailableAcceleration = ""
    GPUbackend = "none"
    Tasmanian.isCudaEnabled() && (GPUbackend = "CUDA")
    Tasmanian.isHipEnabled() && (GPUbackend = "ROCm/HIP")
    Tasmanian.isDpcppEnabled() && (GPUbackend = "oneAPI/DPC++")
    format("                   GPU backend: {1:1s}", GPUbackend)
    for s in ["cpu-blas", "gpu-cublas", "gpu-cuda", "gpu-magma"]
        Tasmanian.isAccelerationAvailable(grid, s) && (AvailableAcceleration *= " " * s)
    end
    format("        Available acceleration:{1:1s}", AvailableAcceleration)

    format("                Available GPUs:")
    if getNumGPUs() > 0
        for GPU in 0:(Tasmanian.getNumGPUs() - 1)
            name = Tasmanian.getGPUName(GPU)
            mem = Tasmanian.getGPUMemory(GPU)
            format("     {1:2d}: {2:20s} with{3:6d}MB RAM", GPU, name, mem)
        end
    else
        format("            none")
    end

    # covers the print() in Julia and C++
    print(grid) # empty grid
    grid = makeGlobalGrid(dimension = 2, outputs = 0, depth = 1, type = "level",
        rule = "gauss-gegenbauer", alpha = 3.0)
    print(grid)
    grid = makeGlobalGrid(dimension = 2, outputs = 0, depth = 1, type = "level",
        rule = "gauss-jacobi", alpha = 3.0, beta = 3.0)
    print(grid)
    grid = makeGlobalGrid(dimension = 2, outputs = 1, depth = 1, type = "level",
        rule = "custom-tabulated", custom_filename = sGaussPattersonTableFile)
    print(grid)
    grid = makeSequenceGrid(
        dimension = 2, outputs = 1, depth = 3, type = "level", rule = "rleja")
    loadExpN2!(grid)
    print(grid)
    grid = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 3)
    setDomainTransform!(grid, [2.0 3.0])
    print(grid)
    grid = makeWaveletGrid(dimension = 2, outputs = 1, depth = 1)
    print(grid)
    grid = makeFourierGrid(dimension = 3, outputs = 1, depth = 1, type = "level")
    enableAcceleration!(grid, "gpu-cuda")
    print(grid)
end

"""
    Test reading and writing of Global grids.
"""
function checkReadWriteGlobal()
    # dimension, outputs, depth, type, rule, alpha, beta, useTransform, loadFunction, limitLevels
    lGrids = [[3, 2, 2, "level", "leja", 0.0, 0.0, false, false, false],
        [2, 1, 4, "level", "clenshaw-curtis", 0.0, 0.0, false, false, false],
        [3, 1, 3, "level", "rleja", 0.0, 0.0, true, false, false],
        [3, 1, 2, "iptotal", "chebyshev", 0.0, 0.0, false, true, false],
        [3, 1, 3, "level", "leja", 0.0, 0.0, true, true, false],
        [3, 1, 3, "level", "leja", 0.0, 0.0, true, true, true],
        [2, 1, 5, "qptotal", "gauss-hermite", 1.0, 3.0, false, false, false],
        [3, 1, 2, "level", "gauss-laguerre", 3.0, 0.0, false, false, true]]
    Transform = [0.0 1.0; 0.0 1.0; -2.0 -1.0]

    for lT in lGrids
        gridA = TasmanianSG()
        gridB = TasmanianSG()

        if lT[8]
            makeGlobalGrid!(
                gridA, dimension = lT[1], outputs = lT[2], depth = lT[3], type = lT[4],
                rule = lT[5], alpha = lT[6], beta = lT[7], level_limits = [3, 2, 1])
        else
            makeGlobalGrid!(gridA, dimension = lT[1], outputs = lT[2], depth = lT[3],
                type = lT[4], rule = lT[5], alpha = lT[6])
        end
        #gridA.print()
        lT[8] && setDomainTransform!(gridA, Transform)
        lT[9] && loadExpN2!(gridA)

        write(gridA, "testSave", binary = false)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        write(gridA, "testSave", binary = true)
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        makeGlobalGrid!(
            gridB, dimension = 1, outputs = 0, depth = 1, type = "level", rule = "rleja")
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        copyGrid!(gridB, gridA)
        compareGrids(gridA, gridB)
    end

    # Test an error message from wrong read.
    try
        read!(gridB, "Test_If_Bogus_Filename_Produces_an_Error")
    catch(e)
        !occursin("Bogus", e.msg) &&
            throws(TasmanianInputError("ERROR in test: Reading a bogus file properly failed, but the error information is wrong."))
    end

    # custom rule test
    gridA = TasmanianSG()
    gridB = TasmanianSG()
    makeGlobalGrid!(gridA, dimension = 2, outputs = 0, depth = 4, type = "level",
        rule = "custom-tabulated", custom_filename = sGaussPattersonTableFile)
    makeGlobalGrid!(gridB, dimension = 2, outputs = 0, depth = 4,
        type = "level", rule = "gauss-patterson")
    compareGrids(gridA, gridB, bTestRuleNames = false)
    write(gridA, "testSave", binary = false)
    makeGlobalGrid!(gridB, dimension = 2, outputs = 0, depth = 4,
        type = "level", rule = "clenshaw-curtis")
    read!(gridB, "testSave")
    compareGrids(gridA, gridB)
    write(gridA, "testSave", binary = true)
    makeGlobalGrid!(
        gridB, dimension = 3, outputs = 0, depth = 4, type = "level", rule = "leja")
    read!(gridB, "testSave")
    compareGrids(gridA, gridB)
end

"""
    Test reading and writing of Sequence grids.
"""
function checkReadWriteSequence()
    # dimension, outputs, depth, type, rule, useTransform, loadFunction, limitLevels
    lGrids = [[3, 2, 2, "level", "leja", false, false, false],
        [2, 1, 4, "level", "max-lebesgue", false, false, false],
        [3, 1, 3, "level", "rleja", true, false, false],
        [3, 1, 3, "level", "rleja", true, false, true],
        [3, 1, 2, "iptotal", "min-delta", false, true, false],
        [3, 1, 3, "level", "leja", true, true, false],
        [3, 1, 3, "level", "leja", true, true, true]]

    for lT in lGrids
        gridA = TasmanianSG()
        gridB = TasmanianSG()

        if lT[8]
            makeSequenceGrid!(gridA, dimension = lT[1], outputs = lT[2], depth = lT[3],
                type = lT[4], rule = lT[5], level_limits = [2, 3, 1])
        else
            makeSequenceGrid!(gridA, dimension = lT[1], outputs = lT[2],
                depth = lT[3], type = lT[4], rule = lT[5])
        end
        lT[6] && setDomainTransform!(gridA, [0.0 1.0; 0.0 1.0; -2.0 -1.0])
        lT[7] && loadExpN2!(gridA)

        write(gridA, "testSave", binary = false)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        write(gridA, "testSave", binary = true)
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        makeGlobalGrid!(
            gridB, dimension = 1, outputs = 0, depth = 1, type = "level", rule = "rleja")
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        copyGrid!(gridB, gridA)
        compareGrids(gridA, gridB)
    end
end

"""
    Test reading and writing of Localp grids.
"""
function checkReadWriteLocalp()
    # dimension, outputs, depth, order, rule, useTransform, loadFunction, limitLevels
    lGrids = [[3, 2, 2, 0, "localp", false, false, false],
        [3, 0, 2, 0, "localp-boundary", false, false, false],
        [2, 1, 4, 1, "semi-localp", false, false, false],
        [3, 1, 3, 2, "localp", true, false, false],
        [3, 1, 2, 3, "localp-zero", false, true, false],
        [3, 1, 2, 3, "localp-zero", false, true, true],
        [3, 1, 3, 4, "semi-localp", true, true, false],
        [3, 1, 3, -1, "semi-localp", true, true, false],
        [3, 1, 3, -1, "semi-localp", true, true, true]]

    for lT in lGrids
        gridA = TasmanianSG()
        gridB = TasmanianSG()

        if lT[7]
            makeLocalPolynomialGrid!(
                gridA, dimension = lT[1], outputs = lT[2], depth = lT[3],
                order = lT[4], rule = lT[5], level_limits = [3, 1, 2])
        else
            makeLocalPolynomialGrid!(gridA, dimension = lT[1], outputs = lT[2],
                depth = lT[3], order = lT[4], rule = lT[5])
        end
        lT[6] && setDomainTransform!(gridA, [0.0 1.0; 0.0 1.0; -2.0 -1.0])
        lT[7] && loadExpN2!(gridA)

        write(gridA, "testSave", binary = false)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        write(gridA, "testSave", binary = true)
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        makeGlobalGrid!(
            gridB, dimension = 1, outputs = 0, depth = 1, type = "level", rule = "rleja")
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        copyGrid!(gridB, gridA)
        compareGrids(gridA, gridB)
    end
end

"""
    Test reading and writing of Wavelet grids.
"""
function checkReadWriteWavelet()
    # dimension, outputs, depth, order, useTransform, loadFunction, level_limits
    lGrids = [Any[3, 2, 2, 1, false, false, false],
        Any[3, 0, 2, 1, false, false, false],
        Any[2, 1, 4, 1, false, false, false],
        Any[3, 1, 1, 3, true, false, false],
        Any[3, 1, 1, 3, true, false, true],
        Any[3, 1, 2, 1, false, true, false],
        Any[3, 1, 2, 3, true, true, true],
        Any[3, 1, 2, 3, true, true, false]]

    for lT in lGrids
        gridA = TasmanianSG()
        gridB = TasmanianSG()

        if lT[7]
            makeWaveletGrid!(gridA, dimension = lT[1], outputs = lT[2],
                depth = lT[3], order = lT[4], level_limits = [1, 1, 2])
        else
            makeWaveletGrid!(
                gridA, dimension = lT[1], outputs = lT[2], depth = lT[3], order = lT[4])
        end

        lT[5] && setDomainTransform!(gridA, [0.0 1.0; 0.0 1.0; -2.0 -1.0])
        lT[6] && loadExpN2!(gridA)

        write(gridA, "testSave", binary = false)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        write(gridA, "testSave", binary = true)
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        makeGlobalGrid!(
            gridB, dimension = 1, outputs = 0, depth = 1, type = "level", rule = "rleja")
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        copyGrid!(gridB, gridA)
        compareGrids(gridA, gridB)
    end
end

"""
    Test reading and writing of Fourier grids.
"""
function checkReadWriteFourier()
    # dimension, outputs, depth, useTransform, loadFunction, useLevelLimits
    lGrids = [Any[3, 2, 2, false, false, false],
        Any[2, 1, 4, false, false, false],
        Any[3, 1, 1, true, false, false],
        Any[3, 1, 1, true, false, true],
        Any[3, 1, 2, false, true, false],
        Any[3, 1, 2, true, true, true],
        Any[3, 1, 2, true, true, false]]

    for lT in lGrids
        gridA = TasmanianSG()
        gridB = TasmanianSG()

        if lT[6]
            makeFourierGrid!(gridA, dimension = lT[1], outputs = lT[2],
                depth = lT[3], type = "level", level_limits = [1, 1, 2])
        else
            makeFourierGrid!(
                gridA, dimension = lT[1], outputs = lT[2], depth = lT[3], type = "level")
        end
        lT[4] && setDomainTransform!(gridA, [0.0 1.0; 0.0 1.0; -2.0 -1.0])
        lT[5] && loadExpN2!(gridA)

        write(gridA, "testSave", binary = false)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        write(gridA, "testSave", binary = true)
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        read!(gridB, "testSave")
        compareGrids(gridA, gridB)

        makeGlobalGrid!(
            gridB, dimension = 1, outputs = 0, depth = 1, type = "level", rule = "rleja")
        makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
        copyGrid!(gridB, gridA)
        compareGrids(gridA, gridB)
    end
end

function checkCopySubgrid_(grids, Grids)
    lValues = [string(g) * "Values" for g in Grids]
    DValues = Dict{String, Vector{Float64}}()
    for v in lValues
        DValues[v] = []
    end

    Points = getPoints(grids[:gridTotal])
    for iI in axes(Points, 2)
        lModel = [i * exp(Points[1, iI] + Points[2, iI]) for i in 1:6]
        DValues["gridTotalValues"] = vcat(DValues["gridTotalValues"], lModel)
        DValues["gridRef1Values"] = vcat(DValues["gridRef1Values"], lModel[1:1])
        DValues["gridRef2Values"] = vcat(DValues["gridRef2Values"], lModel[1:3])
        DValues["gridRef3Values"] = vcat(DValues["gridRef3Values"], lModel[3:5])
        DValues["gridRef4Values"] = vcat(DValues["gridRef4Values"], lModel[6:6])
    end

    for (g, v) in zip(Grids, lValues)
        loadNeededPoints!(grids[Symbol(g)], DValues[v])
    end
    grid = copyGrid(grids[:gridTotal])
    @test grid == grids[:gridTotal]
    grid = copyGrid(grids[:gridTotal], 0, 1)
    @test grid == grids[:gridRef1]
    grid = copyGrid(grids[:gridTotal], 0, 3)
    @test grid == grids[:gridRef2]
    grid = copyGrid(grids[:gridTotal], 2, 5)
    @test grid == grids[:gridRef3]
    grid = copyGrid(grids[:gridTotal], 5, 6)
    @test grid == grids[:gridRef4]
end

function checkCopySubgrid()
    # outputs:   source       0         0, 1, 2     2, 3, 4       5
    Grids = [:gridTotal, :gridRef1, :gridRef2, :gridRef3, :gridRef4]
    Depths = [4, 4, 4, 3, 4]
    Outputs = [6, 1, 3, 3, 1]
    Make = [:makeGlobalGrid,
        :makeSequenceGrid,
        :makeLocalPolynomialGrid,
        :makeWaveletGrid,
        :makeFourierGrid
    ]
    Args = [(depth = 4, type = "level", rule = "clenshaw-curtis"),
        (depth = 4, type = "level", rule = "rleja"),
        (depth = 4, order = 2),
        (depth = 3,),
        (depth = 4, type = "level")
    ]

    grid = Dict()
    for (i, outputs) in enumerate(Outputs)
        grid[Grids[i]] = makeGlobalGrid(dimension = 2, outputs = outputs, depth = 4,
            type = "level", rule = "clenshaw-curtis")
    end
    checkCopySubgrid_(grid, Grids)

    for (i, outputs) in enumerate(Outputs)
        grid[Grids[i]] = makeSequenceGrid(
            dimension = 2, outputs = outputs, depth = 4, type = "level", rule = "rleja")
    end
    checkCopySubgrid_(grid, Grids)

    for (i, outputs) in enumerate(Outputs)
        grid[Grids[i]] = makeLocalPolynomialGrid(
            dimension = 2, outputs = outputs, depth = 4, order = 2)
    end
    checkCopySubgrid_(grid, Grids)

    for (i, outputs) in enumerate(Outputs)
        grid[Grids[i]] = makeWaveletGrid(dimension = 2, outputs = outputs, depth = 3)
    end
    checkCopySubgrid_(grid, Grids)

    for (i, outputs) in enumerate(Outputs)
        grid[Grids[i]] = makeFourierGrid(
            dimension = 2, outputs = outputs, depth = 4, type = "level")
    end
    checkCopySubgrid_(grid, Grids)
end

function checkReadWriteMisc_(gridA)
    gridB = TasmanianSG()
    write(gridA, "testSave", binary = false)
    read!(gridB, "testSave")
    @test compareGrids(gridA, gridB)

    write(gridA, "testSave", binary = true)
    makeLocalPolynomialGrid!(gridB, dimension = 1, outputs = 1, depth = 0)
    read!(gridB, "testSave")
    @test compareGrids(gridA, gridB)

    gridA = makeSequenceGrid(
        dimension = 1, outputs = 1, depth = 0, type = "level", rule = "leja")
    gridB = makeLocalPolynomialGrid(dimension = 1, outputs = 1, depth = 0)
    gridB = copyGrid(gridA)
    @test compareGrids(gridA, gridB)
end

"""
    Test reading and writing of domain transforms and testing all rules.
"""
function checkReadWriteMisc()
    Transform = vcat([0.0 1.0], [0.0 1.0], [-2.0 -1.0])

    gridA = makeGlobalGrid(
        dimension = 3, outputs = 2, depth = 4, type = "level", rule = "clenshaw-curtis")
    setDomainTransform!(gridA, Transform)
    setConformalTransformASIN!(gridA, [3, 4, 5])
    checkReadWriteMisc_(gridA)

    gridA = makeGlobalGrid(
        dimension = 3, outputs = 2, depth = 4, type = "level", rule = "gauss-legendre")
    setConformalTransformASIN!(gridA, [3, 5, 1])
    checkReadWriteMisc_(gridA)

    gridA = makeSequenceGrid(
        dimension = 2, outputs = 2, depth = 5, type = "level", rule = "leja")
    setConformalTransformASIN!(gridA, [0, 4])
    checkReadWriteMisc_(gridA)

    gridA = makeLocalPolynomialGrid(
        dimension = 3, outputs = 1, depth = 4, order = 2, rule = "localp")
    setDomainTransform!(gridA, Transform)
    setConformalTransformASIN!(gridA, [5, 3, 0])
    checkReadWriteMisc_(gridA)

    getNumPoints(gridA)

    # Make a grid with every possible rule (catches false-positive and memory crashes)
    for type in Tasmanian.GlobalTypes
        for rule in Tasmanian.GlobalRules
            if occursin("custom-tabulated", rule)
                gridA = makeGlobalGrid(dimension = 2, outputs = 0, depth = 2, type = type,
                    rule = rule, custom_filename = sGaussPattersonTableFile)
            else
                gridA = makeGlobalGrid(
                    dimension = 2, outputs = 0, depth = 2, type = type, rule = rule)
            end
            write(gridA, "testSave", binary = false)
            gridB = TasmanianSG()
            read!(gridB, "testSave")
            @test compareGrids(gridA, gridB)
            gridB = makeGlobalGrid(dimension = 1, outputs = 0, depth = 0,
                type = "level", rule = "clenshaw-curtis")
            write(gridA, "testSave", binary = true)
            read!(gridB, "testSave")
        end
    end

    for type in Tasmanian.GlobalTypes
        for rule in Tasmanian.SequenceRules
            gridA = makeSequenceGrid(
                dimension = 2, outputs = 1, depth = 3, type = type, rule = rule)
            write(gridA, "testSave", binary = false)
            gridB = TasmanianSG()
            read!(gridB, "testSave")
            @test compareGrids(gridA, gridB)
            gridB = makeGlobalGrid(dimension = 1, outputs = 0, depth = 0,
                type = "level", rule = "clenshaw-curtis")
            write(gridA, "testSave", binary = true)
            read!(gridB, "testSave")
            @test compareGrids(gridA, gridB)
        end
    end
end

"""
    Test reading and writing of the CustomTabulated class.
"""
function checkReadWriteCustomTabulated()
    description = "testCT"
    create_nodes(j) = collect(range(-1.0, 1.0, j))
    function create_weights(j)
        weights = collect(range(0.0, 1.0, j))
        weights = 2.0 * weights / sum(weights)
        return weights
    end
    # Read and write from explicitly given data.
    for i in 0:3
        num_levels = i
        num_nodes = [3 * j for j in 1:i]
        precision = [2 * j - 1 for j in 1:i]
        nodes = [create_nodes(j) for j in num_nodes]
        weights = [create_weights(j) for j in num_nodes]
        ctA = makeCustomTabulatedFromData(
            num_levels, num_nodes, precision, nodes, weights, description)
        ctB = CustomTabulated()
        write(ctA, "testSave")
        read!(ctB, "testSave")
        for j in 1:num_levels
            read_weights, read_nodes = getWeightsNodes(ctA, j - 1)
            @test read_weights ≈ weights[j]
            @test read_nodes ≈ nodes[j]
        end
    end
    # Read and write from a file.
    ctA = makeCustomTabulatedFromFile(sGaussPattersonTableFile)
    ctB = CustomTabulated()
    write(ctA, "testSave")
    read!(ctB, "testSave")
    compareCustomTabulated(ctA, ctB)
    for i in 0:(getNumLevels(ctB) - 1)
        read_weights, read_nodes = getWeightsNodes(ctB, i)
        grid = makeGlobalGrid(
            dimension = 1, outputs = 0, depth = i, type = "level", rule = "gauss-patterson")
        @test read_weights ≈ getQuadratureWeights(grid)
        @test read_nodes' ≈ getPoints(grid)
    end

    # Test an error message from wrong read.
    try
        read!(ctB, "Test_If_Bogus_Filename_Produces_an_Error")
    catch(e)
        !occursin("Bogus", e.msg) &&
            throws(TasmanianInputError("ERROR in test: Reading a bogus file properly failed, but the error information is wrong."))
    end
end

"""
    Test makeGlobalGridCustom(), which creates a grid from a CustomTabulated instance.
"""
function checkGlobalGridCustom()
    gridA = makeGlobalGrid(dimension = 1, outputs = 1, depth = 3, type = "level",
        rule = "custom-tabulated", custom_filename = sGaussPattersonTableFile)
    ct = makeCustomTabulatedFromFile(sGaussPattersonTableFile)
    gridB = makeGlobalGridCustom(
        dimension = 1, outputs = 1, depth = 3, type = "level", ct = ct)
    compareGrids(gridA, gridB)
    gridA = makeGlobalGrid(dimension = 2, outputs = 1, depth = 3, type = "level",
        rule = "custom-tabulated", custom_filename = sGaussPattersonTableFile)
    gridB = makeGlobalGridCustom(
        dimension = 2, outputs = 1, depth = 3, type = "level", ct = ct)
    @test compareGrids(gridA, gridB)
end

@testset verbose=true "Testing core I/O" begin
    @testset "Test the version I/O, available acceleration, gpu info, print stats" checkMetaIO()
    @testset "Test reading and writing global grids" checkReadWriteGlobal()
    @testset "Test reading and writing sequence grids" checkReadWriteSequence()
    @testset "Test reading and writing localp grids" checkReadWriteLocalp()
    @testset "Test reading and writing wavelet grids" checkReadWriteWavelet()
    @testset "Test reading and writing Fourier grids" checkReadWriteFourier()
    @testset "Test grid copy" checkCopySubgrid()
    @testset "Test reading and writing of domain transforms and testing all rules" checkReadWriteMisc()
    @testset "Test reading and writing custom tabulated grids" checkReadWriteCustomTabulated()
    @testset "Test makeGlobalGridCustom(), which creates a grid from a CustomTabulated instance" checkGlobalGridCustom()
end
