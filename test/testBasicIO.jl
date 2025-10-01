using Tasmanian
using Test

include("testCommon.jl")

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
            (depth = 4, type = "level",)
            ]

    for elems in zip(Make, Args)
        make, kwargs = elems
        for (i, output) in enumerate(Outputs)
            grid = Grids[i]
            @eval $grid = $make(;dimension = 2, outputs = $output, $kwargs...)
        end        
        lValues = [string(g)*"Values" for g in Grids]
        DValues = Dict{String, Vector{Float64}}()
        for v in lValues
            DValues[v] = []
        end

        aPoints = getPoints(gridTotal)
        for iI in axes(aPoints, 2)
            lModel = [i * exp(aPoints[1, iI] + aPoints[2, iI]) for i in 1:6]
            DValues["gridTotalValues"] = vcat(DValues["gridTotalValues"], lModel)
            DValues["gridRef1Values"] = vcat(DValues["gridRef1Values"], lModel[1:1])
            DValues["gridRef2Values"] = vcat(DValues["gridRef2Values"], lModel[1:3])
            DValues["gridRef3Values"] = vcat(DValues["gridRef3Values"], lModel[3:5])
            DValues["gridRef4Values"] = vcat(DValues["gridRef4Values"], lModel[6:6])
        end

        for (g, v) in zip(Grids, lValues)
            @eval loadNeededPoints!($g, $(DValues[v]))
        end
        grid = copyGrid(gridTotal)
        @test grid == gridTotal
        grid = copyGrid(gridTotal, 0, 1)
        @test grid == gridRef1
        grid = copyGrid(gridTotal, 0, 3)
        @test grid == gridRef2
        grid = copyGrid(gridTotal, 2, 5)
        @test grid == gridRef3
        grid = copyGrid(gridTotal, 5, 6)
        @test grid == gridRef4
    end
end

#=
    def checkReadWriteMisc(self):
        '''
        Test reading and writing of domain transforms and testing all rules.
        '''
        lGrids = ['gridA.makeGlobalGrid(3, 2, 4, "level", "clenshaw-curtis"); gridA.setDomainTransform(aTransform); gridA.setConformalTransformASIN(np.array([3,4,5]))',
                  'gridA.makeGlobalGrid(3, 2, 4, "level", "gauss-legendre"); gridA.setConformalTransformASIN(np.array([3,5,1]))',
                  'gridA.makeSequenceGrid(2, 2, 5, "level", "leja"); gridA.setConformalTransformASIN(np.array([0,4]))',
                  'gridA.makeLocalPolynomialGrid(3, 1, 4, 2, "localp"); gridA.setDomainTransform(aTransform); gridA.setConformalTransformASIN(np.array([5,3,0]))',
                  'gridA.getNumPoints()']
        aTransform = np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]])

        for sGrid in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            exec(sGrid)
            gridA.write("testSave", bUseBinaryFormat = False)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridB.makeSequenceGrid(1, 1, 0, "level", "leja");
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.copyGrid(gridA)
            ttc.compareGrids(gridA, gridB)

        # Make a grid with every possible rule (catches false-positive and memory crashes)
        for type in TasmanianSG.lsTsgGlobalTypes:
            for rule in TasmanianSG.lsTsgGlobalRules:
                if ("custom-tabulated" in rule):
                    gridA.makeGlobalGrid(2, 0, 2, type, rule, sCustomFilename = tdata.sGaussPattersonTableFile)
                else:
                    gridA.makeGlobalGrid(2, 0, 2, type, rule)
                gridA.write("testSave", bUseBinaryFormat = False)
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)
                gridB.makeGlobalGrid(1, 0, 0, "level", "clenshaw-curtis")
                gridA.write("testSave", bUseBinaryFormat = True)
                gridB.read("testSave")

        for type in TasmanianSG.lsTsgGlobalTypes:
            for rule in TasmanianSG.lsTsgSequenceRules:
                gridA.makeSequenceGrid(2, 1, 3, type, rule)
                gridA.write("testSave", bUseBinaryFormat = False)
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)
                gridB.makeGlobalGrid(1, 0, 0, "level", "clenshaw-curtis")
                gridA.write("testSave", bUseBinaryFormat = True)
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)

    def checkReadWriteCustomTabulated(self):
        '''
        Test reading and writing of the CustomTabulated class.
        '''
        description = "testCT"
        def create_nodes(j):
            return np.linspace(-1.0, 1.0, num=j)
        def create_weights(j):
            weights = np.linspace(0.0, 1.0, num=j)
            weights = 2.0 * weights / weights.sum()
            return weights
        # Read and write from explicitly given data.
        for i in range(4):
            num_levels = i
            num_nodes = np.array([3*(j+1) for j in range(i)])
            precision = np.array([2*(j+1)-1 for j in range(i)])
            nodes = [create_nodes(j) for j in num_nodes]
            weights = [create_weights(j) for j in num_nodes]
            ctA = TasmanianSG.makeCustomTabulatedFromData(num_levels, num_nodes, precision, nodes, weights, description)            # .
            ctB = TasmanianSG.CustomTabulated()
            ctA.write("testSave")
            ctB.read("testSave")
            for i in range(num_levels):
                read_weights, read_nodes = ctB.getWeightsNodes(i)
                np.testing.assert_almost_equal(read_weights, weights[i])
                np.testing.assert_almost_equal(read_nodes, nodes[i])
        # Read and write from a file.
        ctA = TasmanianSG.makeCustomTabulatedFromFile(tdata.sGaussPattersonTableFile)
        ctB = TasmanianSG.CustomTabulated()
        ctA.write("testSave")
        ctB.read("testSave")
        ttc.compareCustomTabulated(ctA, ctB)
        for i in range(ctB.getNumLevels()):
            read_weights, read_nodes = ctB.getWeightsNodes(i)
            grid = TasmanianSG.makeGlobalGrid(1, 0, i, "level", "gauss-patterson")
            np.testing.assert_almost_equal(read_weights, grid.getQuadratureWeights())
            np.testing.assert_almost_equal(read_nodes, grid.getPoints().flatten())
        # Test an error message from wrong read.
        try:
            ctB.read("Test_If_Bogus_Filename_Produces_an_Error")
        except TasmanianSG.TasmanianInputError as TSGError:
            TSGError.bShowOnExit = False
            self.assertEqual(TSGError.sVariable, "sFilename", "Reading a bogus file properly failed, but the error information is wrong.")

    def checkGlobalGridCustom(self):
        '''
        Test makeGlobalGridCustom(), which creates a grid from a CustomTabulated instance.
        '''
        gridA = TasmanianSG.makeGlobalGrid(1, 1, 3, "level", "custom-tabulated", [], 0.0, 0.0, tdata.sGaussPattersonTableFile)
        ct = TasmanianSG.makeCustomTabulatedFromFile(tdata.sGaussPattersonTableFile)
        gridB = TasmanianSG.makeGlobalGridCustom(1, 1, 3, "level", ct)
        ttc.compareGrids(gridA, gridB)
        gridA = TasmanianSG.makeGlobalGrid(2, 1, 3, "level", "custom-tabulated", [], 0.0, 0.0, tdata.sGaussPattersonTableFile)
        gridB = TasmanianSG.makeGlobalGridCustom(2, 1, 3, "level", ct)
        ttc.compareGrids(gridA, gridB)

    def performIOTest(self):
        self.checkMetaIO()

        self.checkReadWriteGlobal()
        self.checkReadWriteSequence()
        self.checkReadWriteLocalp()
        self.checkReadWriteWavelet()
        self.checkReadWriteFourier()

        self.checkCopySubgrid()
        self.checkReadWriteMisc()

        self.checkReadWriteCustomTabulated()
        self.checkGlobalGridCustom()
=#

@testset verbose = true "Testing core I/O test" begin
    @testset "checkCopySubgrid" checkCopySubgrid()
end
