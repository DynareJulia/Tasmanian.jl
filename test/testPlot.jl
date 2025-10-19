using Tasmanian
using Plots
using .RecipesBase
using VisualRegressionTests
using Test

include("../src/plot.jl")

grid = makeLocalPolynomialGrid(
    dimension = 2, outputs = 1, depth = 5, order = 1, rule = "localp")

p1 = scattergridpoints2d(grid)

model(x) = exp.(-x[1, :] .^ 2) .* cos.(x[2, :])
points = getPoints(grid)
values = model(points)
loadNeededValues!(grid, values)
p2 = scattergridresponse2d(grid)

x = -1:0.1:1
y = -1:0.1:1
nx = length(x)
ny = length(y)

points = vcat(vec(ones(ny) * x')',
    vec(y * ones(1, nx))')

p3 = scatterresponse2d(grid, points)

p4 = heatmapresponse2d(grid, x, y)

p5 = scattererrors2d(grid, points, mapslices(model, points, dims = 1))

p6 = heatmaperrors2d(grid, x, y, model)

display(plot(p1, p2, p3, p4, p5, p6, layout = (3, 2)))

@plottest plot(p1, p2, p3, p4, p5, p6, layout = (3, 2)) "refplots/plot1.png"
