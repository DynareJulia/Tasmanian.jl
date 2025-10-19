using Tasmanian
using Plots
using Random

negbox(x) = x * 2 - 1
neg2unit(x) = (x + 1) / 2

function ex1()
    which_basis = 1 #1= linear basis functions -> Check the manual for other options
    tsg = makeLocalPolynomialGrid(
        dimension = 2, outputs = 1, depth = 5, order = which_basis, rule = "localp")

    # sparse grid points from that object
    Points = getPoints(tsg)

    # measure perf at N randomly chosen points
    tfun(x) = cos(0.5 * pi * x[1]) * cos(0.5 * pi * x[2])

    Random.seed!(1)
    N = 1000
    randPnts = negbox.(rand(2, N))

    # truth
    truth = mapslices(tfun, randPnts, dims = 1)

    # values on sparse grid
    spVals = mapslices(tfun, Points, dims = 1)

    # load points needed for such values
    loadNeededPoints!(tsg, spVals)

    # Plots
    p1 = scattergridpoints2d(tsg, m = (:black, 2), title = "test points")
    p2 = scattergridresponse2d(
        tsg, Points[1, :], Points[2, :], m = (:red, 2), title = "sparse grid")
    p3 = surfacegridresponse2d(tsg, randPnts[1, :], randPnts[2, :], title = "prediction")
    p4 = surfaceerrors2d(tsg, randPnts, truth, title = "error")
    p = plot(p1, p2, p3, p4, layout = (2, 2), legend = false)

    display(p)
    return p
end

p = ex1()
