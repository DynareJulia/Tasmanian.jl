@userplot ScatterPoints2D
@userplot ScatterGridPoints2D
@userplot ScatterResponse2D
@userplot ScatterGridResponse2D
@userplot ScatterErrors2D
@userplot SurfaceResponse2D
@userplot SurfaceGridResponse2D
@userplot SurfaceErrors2D
@userplot HeatmapResponse2D
@userplot HeatmapErrors2D

function meshgrid(x, y)
    nx = length(x)
    ny = length(y)
    X = vec(ones(ny) * x')'
    Y = vec(y * ones(1, nx))'
    return vcat(X, Y)
end

@recipe function f(mp::ScatterGridPoints2D)
    grid = mp.args[1]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    gridpoints = getPoints(grid)
    seriestype := :scatter
    markercolor --> :blue
    label := false
    title --> "Grid points"
    return (gridpoints[1, :], gridpoints[2, :])
end

@recipe function f(mp::ScatterPoints2D)
    points = mp.args[1]
    size(points, 1) != 2 && throw(TasmanianInputError("ERROR: points must have 2 rows"))
    seriestype := :scatter
    markercolor --> :blue
    label := false
    title --> "Arbitrary points"
    return (points[1, :], points[2, :])
end

@recipe function f(mp::ScatterGridResponse2D)
    grid = mp.args[1]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    points = getPoints(grid)
    values = getLoadedValues(grid)
    seriestype := :scatter
    markercolor --> :blue
    label := false
    title --> "Grid response"
    return (points[1, :], points[2, :], values[1, :])
end

@recipe function f(mp::ScatterResponse2D)
    grid = mp.args[1]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    points = mp.args[2]
    values = evaluateBatch(grid, points)
    seriestype := :scatter
    markercolor --> :blue
    label := false
    title --> "Approximated function"
    return (points[1, :], points[2, :], values[1, :])
end

@recipe function f(mp::SurfaceGridResponse2D)
    grid = mp.args[1]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    points = getPoints(grid)
    values = getLoadedValues(grid)
    seriestype := :surface
    label := false
    title --> "Grid response"
    return (points[1, :], points[2, :], values[1, :])
end

@recipe function f(mp::SurfaceResponse2D)
    grid = mp.args[1]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    points = mp.args[2]
    values = evaluateBatch(grid, points)
    seriestype := :surface
    label := false
    title --> "Approximated function"
    return (points[1, :], points[2, :], values[1, :])
end

@recipe function f(mp::ScatterErrors2D)
    grid = mp.args[1]
    points = mp.args[2]
    truth = mp.args[3]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    values = evaluateBatch(grid, points)
    errors = abs.(values[1, :]' - truth)
    seriestype := :scatter
    label := false
    title --> "Approximation errors"
    return (points[1, :], points[2, :], errors[1, :])
end

@recipe function f(mp::SurfaceErrors2D)
    grid = mp.args[1]
    points = mp.args[2]
    truth = mp.args[3]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    values = evaluateBatch(grid, points)
    errors = abs.(values[1, :]' - truth)
    seriestype := :surface
    label := false
    title --> "Approximation errors"
    return (points[1, :], points[2, :], errors[1, :])
end

@recipe function f(mp::HeatmapResponse2D)
    grid = mp.args[1]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    x = mp.args[2]
    nx = length(x)
    y = mp.args[3]
    ny = length(y)
    points = vcat([_x for _x in x for _ in y]',
        [_y for _ in x for _y in y]')
    values = evaluateBatch(grid, points)
    mvalues = reshape(values, nx, ny)
    seriestype := :heatmap
    label := false
    title --> "Approximated function"
    return (x, y, mvalues)
end

@recipe function f(mp::HeatmapErrors2D)
    grid = mp.args[1]
    getNumDimensions(grid) != 2 &&
        throw(TasmanianInputError("ERROR: grid must have 2 dimensions"))
    x = mp.args[2]
    nx = length(x)
    y = mp.args[3]
    ny = length(y)
    model = mp.args[4]
    model isa Function ||
        throw(TasmanianInputError("ERROR: the 4th argument must be a function"))
    points = vcat([_x for _x in x for _ in y]',
        [_y for _ in x for _y in y]')
    values = evaluateBatch(grid, points)
    mvalues = reshape(values, nx, ny)
    errors = Matrix{Float64}(undef, nx, ny)
    for i in 1:ny
        for j in 1:nx
            errors[j, i] = abs(mvalues[j, i] - model([x[j], y[i]])[1])
        end
    end
    seriestype := :heatmap
    label := false
    title --> "Approximation error"
    return (x, y, errors)
end

"""
    scatterpoints2d(grid)

if `grid` is a Tasmanian grid, produces a scatter plot of the grid base points.
"""
function scatterpoints2d! end
