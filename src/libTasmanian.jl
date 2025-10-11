using Tasmanian_jll
export Tasmanian_jll

using CEnum

function tsgDestructTasmanianSparseGrid(grid)
    ccall((:tsgDestructTasmanianSparseGrid, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgCopyGrid(destination, source)
    ccall((:tsgCopyGrid, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), destination, source)
end

function tsgCopySubGrid(destination, source, outputs_begin, outputs_end)
    ccall((:tsgCopySubGrid, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cint), destination, source, outputs_begin, outputs_end)
end

# no prototype is found for this function at TasmanianSparseGrid.h:42:13, please use with caution
function tsgGetVersion()
    ccall((:tsgGetVersion, TASlib), Ptr{Cchar}, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:43:13, please use with caution
function tsgGetLicense()
    ccall((:tsgGetLicense, TASlib), Ptr{Cchar}, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:44:5, please use with caution
function tsgGetVersionMajor()
    ccall((:tsgGetVersionMajor, TASlib), Cint, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:45:5, please use with caution
function tsgGetVersionMinor()
    ccall((:tsgGetVersionMinor, TASlib), Cint, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:46:5, please use with caution
function tsgIsOpenMPEnabled()
    ccall((:tsgIsOpenMPEnabled, TASlib), Cint, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:47:5, please use with caution
function tsgIsCudaEnabled()
    ccall((:tsgIsCudaEnabled, TASlib), Cint, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:48:5, please use with caution
function tsgIsHipEnabled()
    ccall((:tsgIsHipEnabled, TASlib), Cint, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:49:5, please use with caution
function tsgIsDpcppEnabled()
    ccall((:tsgIsDpcppEnabled, TASlib), Cint, ())
end

function tsgWrite(grid, filename)
    ccall((:tsgWrite, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}), grid, filename)
end

function tsgWriteBinary(grid, filename)
    ccall((:tsgWriteBinary, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}), grid, filename)
end

function tsgRead(grid, filename)
    ccall((:tsgRead, TASlib), Cint, (Ptr{Cvoid}, Ptr{Cchar}), grid, filename)
end

function tsgMakeGlobalGrid(grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, alpha, beta, custom_filename, limit_levels)
    ccall((:tsgMakeGlobalGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cint}), grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, alpha, beta, custom_filename, limit_levels)
end

function tsgMakeSequenceGrid(grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, limit_levels)
    ccall((:tsgMakeSequenceGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, limit_levels)
end

function tsgMakeLocalPolynomialGrid(grid, dimensions, outputs, depth, order, sRule, limit_levels)
    ccall((:tsgMakeLocalPolynomialGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cint}), grid, dimensions, outputs, depth, order, sRule, limit_levels)
end

function tsgMakeWaveletGrid(grid, dimensions, outputs, depth, order, limit_levels)
    ccall((:tsgMakeWaveletGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cint}), grid, dimensions, outputs, depth, order, limit_levels)
end

function tsgMakeFourierGrid(grid, dimensions, outputs, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgMakeFourierGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, dimensions, outputs, depth, sType, anisotropic_weights, limit_levels)
end

function tsgMakeGridFromCustomTabulated(grid, dimension, outputs, depth, sType, custom_tabulated, anisotropic_weights, limit_levels)
    ccall((:tsgMakeGridFromCustomTabulated, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}), grid, dimension, outputs, depth, sType, custom_tabulated, anisotropic_weights, limit_levels)
end

function tsgUpdateGlobalGrid(grid, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgUpdateGlobalGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, depth, sType, anisotropic_weights, limit_levels)
end

function tsgUpdateSequenceGrid(grid, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgUpdateSequenceGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, depth, sType, anisotropic_weights, limit_levels)
end

function tsgUpdateFourierGrid(grid, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgUpdateFourierGrid, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, depth, sType, anisotropic_weights, limit_levels)
end

function tsgGetAlpha(grid)
    ccall((:tsgGetAlpha, TASlib), Cdouble, (Ptr{Cvoid},), grid)
end

function tsgGetBeta(grid)
    ccall((:tsgGetBeta, TASlib), Cdouble, (Ptr{Cvoid},), grid)
end

function tsgGetOrder(grid)
    ccall((:tsgGetOrder, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumDimensions(grid)
    ccall((:tsgGetNumDimensions, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumOutputs(grid)
    ccall((:tsgGetNumOutputs, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetRule(grid)
    ccall((:tsgGetRule, TASlib), Ptr{Cchar}, (Ptr{Cvoid},), grid)
end

function tsgCopyRuleChars(grid, buffer_size, name, num_actual)
    ccall((:tsgCopyRuleChars, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}), grid, buffer_size, name, num_actual)
end

function tsgGetCustomRuleDescription(grid)
    ccall((:tsgGetCustomRuleDescription, TASlib), Ptr{Cchar}, (Ptr{Cvoid},), grid)
end

function tsgGetNumLoaded(grid)
    ccall((:tsgGetNumLoaded, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumNeeded(grid)
    ccall((:tsgGetNumNeeded, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumPoints(grid)
    ccall((:tsgGetNumPoints, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetLoadedPointsStatic(grid, x)
    ccall((:tsgGetLoadedPointsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgGetLoadedPoints(grid)
    ccall((:tsgGetLoadedPoints, TASlib), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetNeededPointsStatic(grid, x)
    ccall((:tsgGetNeededPointsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgGetNeededPoints(grid)
    ccall((:tsgGetNeededPoints, TASlib), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetPointsStatic(grid, x)
    ccall((:tsgGetPointsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgGetPoints(grid)
    ccall((:tsgGetPoints, TASlib), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetQuadratureWeightsStatic(grid, weights)
    ccall((:tsgGetQuadratureWeightsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, weights)
end

function tsgGetQuadratureWeights(grid)
    ccall((:tsgGetQuadratureWeights, TASlib), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetInterpolationWeightsStatic(grid, x, weights)
    ccall((:tsgGetInterpolationWeightsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, x, weights)
end

function tsgGetInterpolationWeights(grid, x)
    ccall((:tsgGetInterpolationWeights, TASlib), Ptr{Cdouble}, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgLoadNeededPoints(grid, vals)
    ccall((:tsgLoadNeededPoints, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, vals)
end

function tsgLoadNeededValues(grid, vals)
    ccall((:tsgLoadNeededValues, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, vals)
end

function tsgGetLoadedValues(grid)
    ccall((:tsgGetLoadedValues, TASlib), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetLoadedValuesStatic(grid, values)
    ccall((:tsgGetLoadedValuesStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, values)
end

function tsgEvaluate(grid, x, y)
    ccall((:tsgEvaluate, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, x, y)
end

function tsgEvaluateFast(grid, x, y)
    ccall((:tsgEvaluateFast, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, x, y)
end

function tsgIntegrate(grid, q)
    ccall((:tsgIntegrate, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, q)
end

function tsgDifferentiate(grid, x, y)
    ccall((:tsgDifferentiate, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, x, y)
end

function tsgEvaluateBatch(grid, x, num_x, y)
    ccall((:tsgEvaluateBatch, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, num_x, y)
end

function tsgBatchGetInterpolationWeightsStatic(grid, x, num_x, weights)
    ccall((:tsgBatchGetInterpolationWeightsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, num_x, weights)
end

function tsgBatchGetInterpolationWeights(grid, x, num_x)
    ccall((:tsgBatchGetInterpolationWeights, TASlib), Ptr{Cdouble}, (Ptr{Cvoid}, Ptr{Cdouble}, Cint), grid, x, num_x)
end

function tsgIsGlobal(grid)
    ccall((:tsgIsGlobal, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsSequence(grid)
    ccall((:tsgIsSequence, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsLocalPolynomial(grid)
    ccall((:tsgIsLocalPolynomial, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsWavelet(grid)
    ccall((:tsgIsWavelet, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsFourier(grid)
    ccall((:tsgIsFourier, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgSetDomainTransform(grid, a, b)
    ccall((:tsgSetDomainTransform, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, a, b)
end

function tsgIsSetDomainTransfrom(grid)
    ccall((:tsgIsSetDomainTransfrom, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgClearDomainTransform(grid)
    ccall((:tsgClearDomainTransform, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgGetDomainTransform(grid, a, b)
    ccall((:tsgGetDomainTransform, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, a, b)
end

function tsgSetConformalTransformASIN(grid, truncation)
    ccall((:tsgSetConformalTransformASIN, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), grid, truncation)
end

function tsgIsSetConformalTransformASIN(grid)
    ccall((:tsgIsSetConformalTransformASIN, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgClearConformalTransform(grid)
    ccall((:tsgClearConformalTransform, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgGetConformalTransformASIN(grid, truncation)
    ccall((:tsgGetConformalTransformASIN, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), grid, truncation)
end

function tsgClearLevelLimits(grid)
    ccall((:tsgClearLevelLimits, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgGetLevelLimits(grid, limits)
    ccall((:tsgGetLevelLimits, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), grid, limits)
end

function tsgSetAnisotropicRefinement(grid, sType, min_growth, output, level_limits)
    ccall((:tsgSetAnisotropicRefinement, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}), grid, sType, min_growth, output, level_limits)
end

function tsgEstimateAnisotropicCoefficients(grid, sType, output, num_coefficients)
    ccall((:tsgEstimateAnisotropicCoefficients, TASlib), Ptr{Cint}, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}), grid, sType, output, num_coefficients)
end

function tsgEstimateAnisotropicCoefficientsStatic(grid, sType, output, coefficients)
    ccall((:tsgEstimateAnisotropicCoefficientsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}), grid, sType, output, coefficients)
end

function tsgSetGlobalSurplusRefinement(grid, tolerance, output, level_limits)
    ccall((:tsgSetGlobalSurplusRefinement, TASlib), Cvoid, (Ptr{Cvoid}, Cdouble, Cint, Ptr{Cint}), grid, tolerance, output, level_limits)
end

function tsgSetLocalSurplusRefinement(grid, tolerance, sRefinementType, output, level_limits, scale_correction)
    ccall((:tsgSetLocalSurplusRefinement, TASlib), Cvoid, (Ptr{Cvoid}, Cdouble, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}), grid, tolerance, sRefinementType, output, level_limits, scale_correction)
end

function tsgClearRefinement(grid)
    ccall((:tsgClearRefinement, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgMergeRefinement(grid)
    ccall((:tsgMergeRefinement, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgBeginConstruction(grid)
    ccall((:tsgBeginConstruction, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgIsUsingConstruction(grid)
    ccall((:tsgIsUsingConstruction, TASlib), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetCandidateConstructionPointsVoidPntr(grid, sType, output, anisotropic_weights, limit_levels)
    ccall((:tsgGetCandidateConstructionPointsVoidPntr, TASlib), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}), grid, sType, output, anisotropic_weights, limit_levels)
end

function tsgGetCandidateConstructionPointsSurplusVoidPntr(grid, tolerance, sRefType, output, limit_levels, scale_correction)
    ccall((:tsgGetCandidateConstructionPointsSurplusVoidPntr, TASlib), Ptr{Cvoid}, (Ptr{Cvoid}, Cdouble, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}), grid, tolerance, sRefType, output, limit_levels, scale_correction)
end

function tsgGetCandidateConstructionPoints(grid, sType, output, anisotropic_weights, limit_levels, num_points, x)
    ccall((:tsgGetCandidateConstructionPoints, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cdouble}}), grid, sType, output, anisotropic_weights, limit_levels, num_points, x)
end

function tsgGetCandidateConstructionSurplusPoints(grid, tolerance, sRefType, output, limit_levels, scale_correction, num_points, x)
    ccall((:tsgGetCandidateConstructionSurplusPoints, TASlib), Cvoid, (Ptr{Cvoid}, Cdouble, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Ptr{Cdouble}}), grid, tolerance, sRefType, output, limit_levels, scale_correction, num_points, x)
end

function tsgGetCandidateConstructionPointsPythonGetNP(grid, vecx)
    ccall((:tsgGetCandidateConstructionPointsPythonGetNP, TASlib), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), grid, vecx)
end

function tsgGetCandidateConstructionPointsPythonStatic(vecx, x)
    ccall((:tsgGetCandidateConstructionPointsPythonStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), vecx, x)
end

function tsgGetCandidateConstructionPointsPythonDeleteVect(vecx)
    ccall((:tsgGetCandidateConstructionPointsPythonDeleteVect, TASlib), Cvoid, (Ptr{Cvoid},), vecx)
end

function tsgLoadConstructedPoint(grid, x, numx, y)
    ccall((:tsgLoadConstructedPoint, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, numx, y)
end

function tsgFinishConstruction(grid)
    ccall((:tsgFinishConstruction, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgRemovePointsByHierarchicalCoefficient(grid, tolerance, output, scale_correction)
    ccall((:tsgRemovePointsByHierarchicalCoefficient, TASlib), Cvoid, (Ptr{Cvoid}, Cdouble, Cint, Ptr{Cdouble}), grid, tolerance, output, scale_correction)
end

function tsgRemovePointsByHierarchicalCoefficientHardCutoff(grid, num_new, output, scale_correction)
    ccall((:tsgRemovePointsByHierarchicalCoefficientHardCutoff, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}), grid, num_new, output, scale_correction)
end

function tsgEvaluateHierarchicalFunctions(grid, x, num_x, y)
    ccall((:tsgEvaluateHierarchicalFunctions, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, num_x, y)
end

function tsgEvaluateSparseHierarchicalFunctions(grid, x, num_x, pntr, indx, vals)
    ccall((:tsgEvaluateSparseHierarchicalFunctions, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), grid, x, num_x, pntr, indx, vals)
end

function tsgEvaluateSparseHierarchicalFunctionsGetNZ(grid, x, num_x)
    ccall((:tsgEvaluateSparseHierarchicalFunctionsGetNZ, TASlib), Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint), grid, x, num_x)
end

function tsgEvaluateSparseHierarchicalFunctionsStatic(grid, x, num_x, pntr, indx, vals)
    ccall((:tsgEvaluateSparseHierarchicalFunctionsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), grid, x, num_x, pntr, indx, vals)
end

function tsgGetHierarchicalSupportStatic(grid, support)
    ccall((:tsgGetHierarchicalSupportStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, support)
end

function tsgGetHierarchicalCoefficients(grid)
    ccall((:tsgGetHierarchicalCoefficients, TASlib), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetHierarchicalCoefficientsStatic(grid, coeff)
    ccall((:tsgGetHierarchicalCoefficientsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, coeff)
end

function tsgSetHierarchicalCoefficients(grid, c)
    ccall((:tsgSetHierarchicalCoefficients, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, c)
end

function tsgIntegrateHierarchicalFunctions(grid)
    ccall((:tsgIntegrateHierarchicalFunctions, TASlib), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgIntegrateHierarchicalFunctionsStatic(grid, integrals)
    ccall((:tsgIntegrateHierarchicalFunctionsStatic, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, integrals)
end

function tsgPythonGetGlobalPolynomialSpace(grid, interpolation, num_indexes)
    ccall((:tsgPythonGetGlobalPolynomialSpace, TASlib), Ptr{Cint}, (Ptr{Cvoid}, Cint, Ptr{Cint}), grid, interpolation, num_indexes)
end

function tsgGetGlobalPolynomialSpace(grid, interpolation, num_indexes, indexes)
    ccall((:tsgGetGlobalPolynomialSpace, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}), grid, interpolation, num_indexes, indexes)
end

function tsgPrintStats(grid)
    ccall((:tsgPrintStats, TASlib), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgEnableAcceleration(grid, accel)
    ccall((:tsgEnableAcceleration, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}), grid, accel)
end

function tsgEnableAccelerationGPU(grid, accel, gpu)
    ccall((:tsgEnableAccelerationGPU, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint), grid, accel, gpu)
end

function tsgGetAccelerationType(grid)
    ccall((:tsgGetAccelerationType, TASlib), Ptr{Cchar}, (Ptr{Cvoid},), grid)
end

function tsgSetGPUID(grid, gpuID)
    ccall((:tsgSetGPUID, TASlib), Cvoid, (Ptr{Cvoid}, Cint), grid, gpuID)
end

function tsgGetGPUID(grid)
    ccall((:tsgGetGPUID, TASlib), Cint, (Ptr{Cvoid},), grid)
end

# no prototype is found for this function at TasmanianSparseGrid.h:150:5, please use with caution
function tsgGetNumGPUs()
    ccall((:tsgGetNumGPUs, TASlib), Cint, ())
end

function tsgGetGPUMemory(gpu)
    ccall((:tsgGetGPUMemory, TASlib), Cint, (Cint,), gpu)
end

function tsgIsAccelerationAvailable(accel)
    ccall((:tsgIsAccelerationAvailable, TASlib), Cint, (Ptr{Cchar},), accel)
end

function tsgGetGPUName(gpu, num_buffer, buffer, num_actual)
    ccall((:tsgGetGPUName, TASlib), Cvoid, (Cint, Cint, Ptr{Cchar}, Ptr{Cint}), gpu, num_buffer, buffer, num_actual)
end

function tsgDeleteInts(p)
    ccall((:tsgDeleteInts, TASlib), Cvoid, (Ptr{Cint},), p)
end

# no prototype is found for this function at TasmanianSparseGrid.h:155:7, please use with caution
function tsgConstructCustomTabulated()
    ccall((:tsgConstructCustomTabulated, TASlib), Ptr{Cvoid}, ())
end

function tsgDestructCustomTabulated(ct)
    ccall((:tsgDestructCustomTabulated, TASlib), Cvoid, (Ptr{Cvoid},), ct)
end

function tsgWriteCustomTabulated(ct, filename)
    ccall((:tsgWriteCustomTabulated, TASlib), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}), ct, filename)
end

function tsgReadCustomTabulated(ct, filename)
    ccall((:tsgReadCustomTabulated, TASlib), Cint, (Ptr{Cvoid}, Ptr{Cchar}), ct, filename)
end

function tsgGetNumLevelsCustomTabulated(ct)
    ccall((:tsgGetNumLevelsCustomTabulated, TASlib), Cint, (Ptr{Cvoid},), ct)
end

function tsgGetNumPointsCustomTabulated(ct, level)
    ccall((:tsgGetNumPointsCustomTabulated, TASlib), Cint, (Ptr{Cvoid}, Cint), ct, level)
end

function tsgGetIExactCustomTabulated(ct, level)
    ccall((:tsgGetIExactCustomTabulated, TASlib), Cint, (Ptr{Cvoid}, Cint), ct, level)
end

function tsgGetQExactCustomTabulated(ct, level)
    ccall((:tsgGetQExactCustomTabulated, TASlib), Cint, (Ptr{Cvoid}, Cint), ct, level)
end

function tsgGetDescriptionCustomTabulated(ct)
    ccall((:tsgGetDescriptionCustomTabulated, TASlib), Ptr{Cchar}, (Ptr{Cvoid},), ct)
end

function tsgGetWeightsNodesStaticCustomTabulated(ct, level, w, x)
    ccall((:tsgGetWeightsNodesStaticCustomTabulated, TASlib), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), ct, level, w, x)
end

function tsgMakeCustomTabulatedFromData(cnum_levels, cnum_nodes, cprecision, cnodes, cweights, cdescription)
    ccall((:tsgMakeCustomTabulatedFromData, TASlib), Ptr{Cvoid}, (Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}), cnum_levels, cnum_nodes, cprecision, cnodes, cweights, cdescription)
end

function tsgGetSubrules(ct, start_index, stride, description)
    ccall((:tsgGetSubrules, TASlib), Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Cint, Ptr{Cchar}), ct, start_index, stride, description)
end

