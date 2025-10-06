module Tasmanian

import Base: show, read, run, write
using LinearAlgebra
using RecipesBase
using Random
using SparseArrays
using Tasmanian_jll

const global TASlib = libtasmaniansparsegrid

# includes
include("libTasmanian.jl")
include("plot.jl")
include("TSG.jl")
include("../examples/examples.jl")

export TasmanianSG, beginConstruction, clearConformalTransform!, clearDomainTransform!, clearLevelLimits!, clearRefinement!, compareGrids, copyGrid,
    differentiate, differentiate!, enableAcceleration!,
    estimateAnisotropicCoefficients, evaluate, evaluateBatch,
    evaluateBatch!, evaluateHierarchicalFunctions, evaluateSparseHierarchicalFunctions,
    evaluateThreadSafe, finishConstruction!, getAlpha, getAccelerationType, getAnisotropicRefinement!, getBeta, getCandidateConstructionPoints,
    getCandidateConstructionPointsSurplus, getConformalTransformASIN, getCustomRuleDescription, getDims,
    getDomainTransform, getHierarchicalCoefficients,
    getInterpolationWeights, getInterpolationWeightsBatch, getGPUID, getGPUMemory, getHierarchicalSupport, getLevelLimits, getLoadedPoints, getLoadedValues,
    getNeededPoints, getNout, getNumDimensions, getGPUName, getNumGPUs, getNumLoaded,
    getNumNeeded, getNumOutputs, getNumPoints, getOrder, getPoints,
    getQuadratureWeights, getRule, getSurplusRefinement!, integrate, integrateHierarchicalFunctions, isAccelerationAvailable, isFourier, isGlobal,
    isLocalPolynomial, isSequence, isSetDomainTransform, isSetConformalTransformASIN, isUsingConstruction, 
    isWavelet, loadConstructedPoint!,
    loadNeededPoints!, loadNeededValues!, makeFourierGrid, makeFourierGrid!, makeGlobalGrid, makeGlobalGrid!, makeLocalPolynomialGrid,
    makeLocalPolynomialGrid!, makeSequenceGrid, makeSequenceGrid!, makeWaveletGrid, makeWaveletGrid!, mergeRefinement!,
    removePointsByHierarchicalCoefficient!, setAnisotropicRefinement!, setConformalTransformASIN!, makeWaveletGrid!,
    setDomainTransform!, setGPUID!, setHierarchicalCoefficients!, setSurplusRefinement!, updateFourierGrid!, updateGlobalGrid!,
    updateLocalPolynomialGrid!, updateSequenceGrid!, updateWaveletGrid!


end # module
