using Tasmanian
using Plots
using Random

negbox(x) = x * 2 - 1
neg2unit(x) = (x + 1) / 2

function ex2(; save_gif = false)
    tol = 1e-5
    K = 7  # max refinement steps
    which_basis = 1 #1= linear basis functions -> Check the manual for other options
    tsg = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 2, order = which_basis, rule = "localp")

    # sparse grid points from that object
    spPoints = getPoints(tsg)

    # test fun 
    tfun(x) = exp(-x[1]^2) * cos(x[2])

    Random.seed!(2)
    N = 1000
    randPnts = negbox.(rand(2, N))
    # truth
    truth = mapslices(tfun, randPnts, dims = 1)

    # values on sparse grid
    spVals = mapslices(tfun, spPoints, dims = 1)
    # load points needed for such values
    loadNeededPoints!(tsg, spVals)

    # evaluate interpolation
    res = evaluateBatch(tsg, randPnts)

    numpoints = size(spPoints, 2)

    @info("error on initial grid:    $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")

    # refinefment loop
    anim = @animate for k in 1:K
        setSurplusRefinement!(tsg, tolerance = tol, refinement_type = "classic")
        if getNumNeeded(tsg) > 0
            spPoints = getNeededPoints(tsg)   # additional set of points required after refinement
            spVals = mapslices(tfun, spPoints, dims = 1)
            # load points needed for such values
            loadNeededPoints!(tsg, spVals)
            numpoints = +size(spPoints, 2)

            # evaluate interpolation
            res = evaluateBatch(tsg, randPnts)
            pred = evaluateBatch(tsg, spPoints)  # prediction on spGrid
            @info("refinement level $k error: $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")

            # plot
            zerone = (-1.1, 1.1)
            p1 = scattergridpoints2d(
                tsg, title = "level $k grid:\n $numpoints points", m = (:black, 1, :+),
                aspect_ratio = :equal, xlims = zerone, ylims = zerone)
            p2 = scatterresponse2d(tsg, randPnts,
                title = "max error: $(round(maximum(abs,res .- truth),digits = 5))",
                m = (:red, 1),
                xlims = zerone, ylims = zerone, zlims = (0, 1))
            p3 = scattergridresponse2d(tsg, title = "grid prediction", m = (:red, 1, 0.2),
                zlims = (0, 1), xlims = zerone, ylims = zerone, zgrid = :black)
            p = plot(p1, p2, p3, layout = (1, 3), leg = false)
        end
    end
    if save_gif
        gif(anim, joinpath(dirname(@__FILE__), "ex2.gif"), fps = 1)
    end
    return p
end

p = ex2()
