using JET
using Format, RecipesPipeline, Unzip

@report_call ignored_modules = (Base, Format) checkMeta()
@report_call ignored_modules = (Base, Format) checkMultiGPU()
@report_call ignored_modules = (Base, Format) checkEvaluateConsistency()

@report_call ignored_modules = (Base, Format) checkMetaIO()
@report_call ignored_modules = (Base, Format) checkReadWriteGlobal()
@report_call ignored_modules = (Base, Format) checkReadWriteSequence()
@report_call ignored_modules = (Base, Format) checkReadWriteLocalp()
@report_call ignored_modules = (Base, Format) checkReadWriteWavelet()
@report_call ignored_modules = (Base, Format) checkReadWriteFourier()
@report_call ignored_modules = (Base, Format) checkCopySubgrid()
@report_call ignored_modules = (Base, Format) checkReadWriteMisc()
@report_call ignored_modules = (Base, Format) checkReadWriteCustomTabulated()
@report_call ignored_modules = (Base, Format) checkGlobalGridCustom()

@report_call ignored_modules = (Base, Format) getSparseGridTests()

@report_call ignored_modules = (Base, Format) checkPathsVersions()
@report_call ignored_modules = (Base, Format) checkMakeAgainstKnown()
@report_call ignored_modules = (Base, Format) checkUpdate()
@report_call ignored_modules = (Base, Format) checkDefaults()
@report_call ignored_modules = (Base, Format) checkLevelLimits()

@report_call ignored_modules = (Base, Format) checkDerivatives()

@report_call ignored_modules = (Base, Format) checkSetClear()      
@report_call ignored_modules = (Base, Format) checkAnisoCoeff()    
@report_call ignored_modules = (Base, Format) checkLocalpSurplus() 
@report_call ignored_modules = (Base, Format) checkFileIO()        
@report_call ignored_modules = (Base, Format) checkConstruction()  
@report_call ignored_modules = (Base, Format) checkRemovePoints()  

@report_call ignored_modules = (Base, Format) checkAgainsKnown()
@report_call ignored_modules = (Base, Format) checkSetCoeffsMergeRegine()
@report_call ignored_modules = (Base, Format) checkIntegrals()
@report_call ignored_modules = (Base, Format) checkValues()

@report_call ignored_modules = (Base, Format, Plots, PlotUtils, RecipesBase,
                                RecipesPipeline, Unzip) ex1()

@report_call ignored_modules = (Base, Format, Plots, PlotUtils, RecipesBase,
                                RecipesPipeline, Unzip) ex2()
@report_call ignored_modules = (Base, Format, Plots, PlotUtils, RecipesBase,
                                RecipesPipeline, Unzip) ex3()
@report_call ignored_modules = (Base, Format) example_01()
@report_call ignored_modules = (Base, Format) example_02()
@report_call ignored_modules = (Base, Format) example_03()
@report_call ignored_modules = (Base, Format) example_04()
@report_call ignored_modules = (Base, Format) example_05()


report_package("Tasmanian", ignored_modules = (Base, Format, Plots, PlotUtils, RecipesBase,
                                RecipesPipeline, Unzip))
