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

using Format
using Tasmanian

function example_03()
    println("\n---------------------------------------------------------------------------------------------------\n")
    println("Example 3: integrate exp(-x1^2 - x2^2) * cos(x3) * cos(x4)")
    println("           for x1, x2 in [-5,5]; x3, x4 in [-2,3]")
    println("           using different rules and total degree polynomial space\n")

    function make_grid(iPrecision, sRule)
        grid = makeGlobalGrid(
            dimension = 4, outputs = 0, depth = iPrecision, type = "qptotal", rule = sRule)
        setDomainTransform!(grid, vcat([-5.0 5.0], [-5.0 5.0], [-2.0 3.0], [-2.0 3.0]))
        return grid
    end

    function print_error(grid)
        fExactIntegral = 1.861816427518323e+00 * 1.861816427518323e+00
        aPoints = getPoints(grid)
        aWeights = getQuadratureWeights(grid)

        fApproximateIntegral = sum(aWeights .*
                                   exp.(-aPoints[1, :] .^ 2 - aPoints[2, :] .^ 2)
                                   .*
                                   cos.(aPoints[3, :]) .* cos.(aPoints[4, :]))
        fError = abs.(fApproximateIntegral - fExactIntegral)
        return format("{1:>10d}{2:>10.2e}", getNumPoints(grid), fError)
    end

    println("               Clenshaw-Curtis      Gauss-Legendre    Gauss-Patterson")
    println(" precision    points     error    points     error    points    error")

    for prec in range(5, 40, step = 5)
        println(format("{1:>10d}{2:1s}{3:1s}{4:1s}",
            prec,
            print_error(make_grid(prec, "clenshaw-curtis")),
            print_error(make_grid(prec, "gauss-legendre-odd")),
            print_error(make_grid(prec, "gauss-patterson"))))
    end

    println("\nAt 311K points the Gauss-Legendre error is O(1.E-1),")
    println("                   Clenshaw-Curtis error is O(1.E-7) at 320K points.")
    println("At 70K points the Gauss-Patterson error is O(1.E-4),")
    println("                  Clenshaw-Curtis needs 158K points to achieve the same.")
end

example_03()
