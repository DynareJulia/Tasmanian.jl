##############################################################################################################################################################################
# Copyright (c) 2017, Miroslav Stoyanov
#
# This file is part of
# Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
#    and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
# THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
# COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
# THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
# IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
##############################################################################################################################################################################

# Need to use parallel version of loadNeededValues provided in AddOn

using Format
using Tasmanian

function  example_05()

    println("\n---------------------------------------------------------------------------------------------------\n")
    println("Example 5: interpolate f(x,y) = exp(-x^2) * cos(y), using leja rule")
    println("           employ adaptive refinement to increase accuracy per samples")

    iNumInputs = 2
    iNumOutputs = 1
    model(aX) = exp(-aX[1] * aX[1]) * cos(aX[2])

    iTestGridSize = 33
    dx = range(-1.0, 1.0, iTestGridSize) # sample on a uniform grid

    aMeshX = [x for y in dx, x in dx]
    aMeshY = [y for y in dx, x in dx]  

    aTestPoints = vcat(vec(aMeshX)', vec(aMeshY)')

    aReferenceValues = exp.(-aTestPoints[1,:].^2) .*  cos.(aTestPoints[2,:])

    function testGrid(grid, aTestPoints, aReferenceValues)
        aResult = evaluateBatch(grid, aTestPoints)
        return maximum(abs.(aResult' - aReferenceValues))
    end

    iInitialLevel = 5

    grid_isotropic = makeGlobalGrid(dimension = iNumInputs, outputs = iNumOutputs, depth = iInitialLevel, type = "level", rule = "leja")
    println(grid_isotropic)
    grid_iptotal = copyGrid(grid_isotropic)
    grid_icurved = copyGrid(grid_isotropic)
    grid_surplus = copyGrid(grid_isotropic)

    iNumThreads = 1
    iBudget = 100
    println(format("{1:>22s}{2:>22s}{3:>22s}{4:>22s}", "isotropic", "iptotal", "ipcurved", "surplus"))
    println(format("{1:>8s}{2:>14s}{1:>8s}{2:>14s}{1:>8s}{2:>14s}{1:>8s}{2:>14s}", "points", "error"))

    bBelowBudget = true
    while(bBelowBudget)
        sInfo = ""
        if getNumLoaded(grid_isotropic) < iBudget
            loadNeededValues!(grid_isotropic, mapslices(model, getNeededPoints(grid_isotropic), dims = 1))
            sInfo *= format("{1:>8d}{2:>14.4e}", getNumLoaded(grid_isotropic),
                            testGrid(grid_isotropic, aTestPoints, aReferenceValues))

            iLevel = 0
            while(getNumNeeded(grid_isotropic) == 0)
                updateGlobalGrid!(grid_isotropic, depth = iLevel, type = "level")
                iLevel += 1
            end
        else
            sInfo *= format("{1:>22s}", "")
        end
        
        if getNumLoaded(grid_iptotal) < iBudget
            loadNeededValues!(grid_iptotal, mapslices(model, getNeededPoints(grid_iptotal), dims = 1))
            sInfo *= format("{1:>8d}{2:>14.4e}", getNumLoaded(grid_iptotal),
                            testGrid(grid_iptotal, aTestPoints, aReferenceValues))

            setAnisotropicRefinement!(grid_iptotal, type = "iptotal", min_growth = 10, output = 0)
        else
            sInfo *= format("{1:>22s}", "")
        end
        
        if getNumLoaded(grid_icurved) < iBudget
            loadNeededValues!(grid_icurved, mapslices(model, getNeededPoints(grid_icurved), dims = 1))
            sInfo *= format("{1:>8d}{2:>14.4e}", getNumLoaded(grid_icurved),
                            testGrid(grid_icurved, aTestPoints, aReferenceValues))

            setAnisotropicRefinement!(grid_icurved, type = "ipcurved", min_growth = 10, output = 0)
        else
            sInfo *= format("{1:>22s}", "")
        end
            
        if getNumLoaded(grid_surplus) < iBudget
            loadNeededValues!(grid_surplus, mapslices(model, getNeededPoints(grid_surplus), dims = 1))
            sInfo *= format("{1:>8d}{2:>14.4e}", getNumLoaded(grid_surplus),
                            testGrid(grid_surplus, aTestPoints, aReferenceValues))

            setSurplusRefinement!(grid_surplus, tolerance = 1.E-8, output =  0)
        else
            sInfo *= format("{1:>22s}", "")
        end
            
        println(sInfo)
        bBelowBudget = (getNumLoaded(grid_isotropic) < iBudget
                       || getNumLoaded(grid_icurved) < iBudget
                       || getNumLoaded(grid_icurved) < iBudget
                        || getNumLoaded(grid_surplus) < iBudget)
    end
end
    

example_05()
