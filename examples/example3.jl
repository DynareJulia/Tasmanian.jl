using Tasmanian
using Plots
using Random

negbox(x) = x * 2 - 1
neg2unit(x) = (x + 1) / 2

function ex3(; save_gif = false)
    tfun(x) = 1.0 / (abs(0.5 - x[1]^4 - x[2]^4) + 0.1)
    tol = 1e-5
    K = 9  # max refinement steps

    which_basis = 1 #1= linear basis functions -> Check the manual for other options
    tsg = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 2, order = which_basis, rule = "localp")

    # domain is [0,1] here
    setDomainTransform!(tsg, [0 1.0; 0 1])

    # sparse grid points from that object
    spPoints = getPoints(tsg)
    Random.seed!(2)
    N = 1000
    randPnts = rand(2, N)
    # truth
    truth = mapslices(tfun, randPnts, dims = 1)

    # values on sparse grid
    spVals = mapslices(tfun, spPoints, dims = 1)
    # load points needed for such values
    loadNeededPoints!(tsg, spVals)

    # evaluate interpolation
    res = evaluateBatch(tsg, randPnts)

    numpoints = size(spPoints, 2)

    @info("error on initial grid:    $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")

    # refinefment loop
    anim = @animate for k in 1:K
        setSurplusRefinement!(tsg, tolerance = tol, refinement_type = "classic")
        if getNumNeeded(tsg) > 0
            spPoints = getNeededPoints(tsg)
            spVals = mapslices(tfun, spPoints, dims = 1)
            # load points needed for such values
            loadNeededPoints!(tsg, spVals)
            numpoints = +size(spPoints, 2)

            # evaluate interpolation
            res = evaluateBatch(tsg, randPnts)
            pred = evaluateBatch(tsg, spPoints)  # prediction on spGrid
            zerone = (-0.1, 1.1)
            @info("refinement level $k error: $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")
            p1 = scattergridpoints2d(
                tsg, title = "level $k grid: $numpoints points", m = (:black, 2, :+),
                aspect_ratio = :equal, xlims = zerone, ylims = zerone)
            p2 = scatterresponse2d(tsg, randPnts, res,
                title = "max error: $(round(maximum(abs,res .- truth),digits = 5))",
                m = (:black, 1))#, zlims=(0,10), camera=(30,70), xlims=zerone, ylims=zerone),
            p3 = scattergridresponse2d(
                tsg, spPoints, title = "grid prediction", m = (:red, 1), zlims = (0, 10),
                camera = (30, 70), xlims = zerone, ylims = zerone, zforeground_color_grid = :black)
            p = plot(p1, p2, p3, layout = (1, 3), leg = false)
            display(p)
        end
        if save_gif
            gif(anim, joinpath(dirname(@__FILE__), "ex3.gif"), fps = 1)
        end
    end
    return p
end

ex3()
