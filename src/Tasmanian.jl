module Tasmanian

import Base: show, read, read!, run, write
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

export TasmanianSG, CustomTabulated, beginConstruction!,
    clearConformalTransform!, clearDomainTransform!,
    clearLevelLimits!, clearRefinement!, compareGrids, copyGrid, copyGrid!,
    differentiate, differentiate!, enableAcceleration!,
    estimateAnisotropicCoefficients, evaluate, evaluateBatch,
    evaluateBatch!, evaluateHierarchicalFunctions,
    evaluateSparseHierarchicalFunctions, evaluateThreadSafe,
    finishConstruction!, getAlpha, getAccelerationType,
    getAnisotropicRefinement!, getBeta,
    getCandidateConstructionPoints,
    getCandidateConstructionPointsSurplus, getConformalTransformASIN,
    getCustomRuleDescription, getDescription, getDims,
    getDomainTransform, getHierarchicalCoefficients,
    getInterpolationWeights, getInterpolationWeightsBatch, getGPUID,
    getGPUMemory, getHierarchicalSupport, getIExact, getLevelLimits, getLicense,
    getLoadedPoints, getLoadedValues, getNeededPoints, getNout,
    getNumDimensions, getNumLevels, getGPUName, getNumGPUs,
    getNumLoaded, getNumNeeded, getNumOutputs, getNumPoints, getOrder,
    getPoints, getQuadratureWeights, getQExact, getRule,
    getSurplusRefinement!, getVersion, getVersionMajor, getVersionMinor, getWeightsNodes, integrate,
    integrateHierarchicalFunctions, isAccelerationAvailable,
    isFourier, isGlobal, isLocalPolynomial, isSequence,
    isSetDomainTransform, isSetConformalTransformASIN,
    isUsingConstruction, isWavelet, loadConstructedPoint!,
    loadNeededPoints!, loadNeededValues!,
    makeCustomTabulatedFromData, makeCustomTabulatedFromFile,
    makeCustomTabulatedSubset, makeFourierGrid, makeFourierGrid!,
    makeGlobalGrid, makeGlobalGrid!, makeGlobalGridCustom, makGlobalGridCustom!, makeLocalPolynomialGrid,
    makeLocalPolynomialGrid!, makeSequenceGrid, makeSequenceGrid!,
    makeWaveletGrid, makeWaveletGrid!, mergeRefinement!, makeWaveletGrid!, read!,
    removePointsByHierarchicalCoefficient!, setAnisotropicRefinement!,
    setConformalTransformASIN!, 
    setDomainTransform!, setGPUID!, setHierarchicalCoefficients!,
    setSurplusRefinement!, updateFourierGrid!, updateGlobalGrid!,
    updateLocalPolynomialGrid!, updateSequenceGrid!,
    updateWaveletGrid!, write


end # module
