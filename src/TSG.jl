GlobalRules = ["clenshaw-curtis", "clenshaw-curtis-zero", "chebyshev", "chebyshev-odd", "gauss-legendre", "gauss-legendre-odd", "gauss-patterson", "leja", "leja-odd", "rleja", "rleja-odd", "rleja-double2", "rleja-double4", "rleja-shifted", "rleja-shifted-even", "rleja-shifted-double", "max-lebesgue", "max-lebesgue-odd", "min-lebesgue", "min-lebesgue-odd", "min-delta", "min-delta-odd", "gauss-chebyshev1", "gauss-chebyshev1-odd", "gauss-chebyshev2", "gauss-chebyshev2-odd", "fejer2", "gauss-gegenbauer", "gauss-gegenbauer-odd", "gauss-jacobi", "gauss-jacobi-odd", "gauss-laguerre", "gauss-laguerre-odd", "gauss-hermite", "gauss-hermite-odd", "custom-tabulated"]
GlobalTypes = ["level", "curved", "iptotal", "ipcurved", "qptotal", "qpcurved", "hyperbolic", "iphyperbolic", "qphyperbolic", "tensor", "iptensor", "qptensor"]
CurvedTypes = ["curved", "ipcurved", "qpcurved"]
RefineTypes = ["classic", "parents", "direction", "fds", "stable"]
SequenceRules = ["leja", "rleja", "rleja-shifted", "max-lebesgue", "min-lebesgue", "min-delta"]
LocalRules = ["localp", "semi-localp", "localp-zero", "localp-boundary"]
AccelTypes = ["none", "cpu-blas", "gpu-default", "gpu-cublas", "gpu-cuda", "gpu-rocblas", "gpu-hip", "gpu-magma"]

struct TasmanianInputError <: Exception
    msg::String
end

Base.showerror(io::IO, e::TasmanianInputError) = print(io, e.msg)
               
function show(io::IO,TSG::TasmanianSG)
    ccall((:tsgPrintStats,TASlib),Nothing,(Ptr{Nothing},),TSG.pGrid)
end

function del(TSG::TasmanianSG)
  ccall((:tsgDestructTasmanianSparseGrid,TASlib),Nothing,(Ptr{Nothing},),TSG.pGrid)
end

"""    getVersion()

returns Tasmanian library hardcoded version
"""
getVersion() = VersionNumber(split(unsafe_string(tsgGetVersion()))[1])

"""
    getLicense()

returns Tasmanian library hardcoded license string
"""
getLicense() = unsafe_string(tsgGetLicense())

"""
    get_VersionMajor()

returns the hardcoded version major int
"""
get_VersionMajor()  = tsgGetVersionMajor()

"""
    get_VersionMinor()

returns the hardcoded version minor int
"""
get_VersionMinor()  = tsgGetVersionMinor()

"""
    isOpenMPEnabled()

returns `true` if the library has been computed with OpenMP support
"""
isOpenMPEnabled() = convert(Bool, tsgIsOpenMPEnabled())

"""
    isCudaEnabled()

returns `true` if the library has been computed with Cuda support
"""
isCudaEnabled() = convert(Bool, tsgIsCudaEnabled())

"""
    isHipEnabled()

returns `true` if the library has been computed with HIP support
"""
isHipEnabled() = convert(Bool, tsgIsHipEnabled())

"""
    isDpcppEnabled()

returns `true` if the library has been computed with DPC++ support
"""
isDpcppEnabled() = convert(tsgIsDpcppEnabled())

"""
    read(tsg::TasmanianSG, filename)

reads the `tsg` grid from a file
discards any existing grid held by `tsg`

`tsg`: an existing grid
`filename`: string indicating a grid file where a grid was
            already written using write from Julia or any other
            Tasmanian interface

output: Bool
        `true`: the read was successful
        `false`: the read failed,
                 check for an error message
"""
function read(tsg::TasmanianSG, filename)
    res = convert(Bool, tsgRead(tsg.pGrid, filename))
    if !res
        throw(TasmanianInputError("ERROR: $filename does not appeaar to be a valid Tasmanian file."))
    end
    return res
end

"""
    write(tsg::TasmanianSG, filename; binary::Bool=true)

writes the `tsg` grid to a file

`tsg`: an existing grid
`filename`: string indicating a grid file where a grid will be written
`binary`: Bool
    `true`: write to a binary file
    `false`: write to an ASCII file
"""
function write(tsg::TasmanianSG, filename; binary::Bool=true)
    if binary
        tsgWriteBinary(tsg.pGrid, filename)
    else
        tsgWrite(tsg.pGrid, filename)
    end
    return nothing
end

function check_anisotropic_weights_(anisotropic_weights, dimension::Int, type::String)
    nweights = (type in CurvedTypes) ? 2*dimension : dimension
    if !isempty(anisotropic_weights)
        if length(anisotropic_weights) != nweights
            throw(TasmanianInputError("ERROR: wrong number of anisotropic_weights, type `$type` needs $nweights weights but length(anisotropic_weights) == $(length(anisotropic_weights))`"))
        end
        if eltype(anisotropic_weights) != Int32
            return convert(Vector{Int32}, anisotropic_weights)
        else
            return anisotropic_weights
        end
    else
        return C_NULL
    end
end

function check_level_limits_(level_limits, dimension)
    n = length(level_limits)
    if n > 0
        if n != dimension
            throw(TasmanianInputError("ERROR: invalid number of level_limits. level_limits needs to have $dimension elements"))
        end
        if eltype(level_limits) != Int32 
            return convert(Vector{Int32}, level_limits)
        else
            return level_limits
        end
    else
        return C_NULL
    end
end

"""
    makeGlobalGrid(; dimensions, outputs, depth, type, rule, anisotropic_weights=[], alpha=0.0, beta=0.0, custom_filename="", level_limits=[])

creates a new sparse grid using a global polynomial rule

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

type: string identifying the tensor selection strategy
     `level`     `curved`     `hyperbolic`     `tensor`
     `iptotal`   `ipcurved`   `iphyperbolic`   `iptensor`
     `qptotal`   `qpcurved`   `qphyperbolic`   `qptensor`

rule: string (defines the 1-D rule that induces the grid)

       Interpolation rules

          Note: the quadrature induced by those rules is constructed
                by integrating the interpolant

          `clenshaw-curtis`    `clenshaw-curtis-zero`      `fejer2`
          `rleja`    `rleja-odd`  `rleja-double2`   `rleja-double4`
          `rleja-shifted`   `rleja-shifted-even`
          `max-lebesgue`    `max-lebesgue-odd`
          `min-lebesgue`    `min-lebesgue-odd`
          `leja`            `leja-odd`
          `min-delta`       `min-delta-odd`

          `chebyshev`       `chebyshev-odd`
            approximation using roots of Chebyshev polynomials
            non-nested case (in contrast to Clenshaw-Curtis nodes)

       Quadrature rules, the weights target exactness with respect
                         to the highest polynomial degree possible

           `gauss-legendre`  `gauss-legendre-odd`
            approximation using roots of polynomials orthogonal in
            measure Uniform

           `gauss-patterson`  (a.k.a. nested Gauss-Legendre)
            Note: the nodes and weights are hard-coded hence there
            is a limit on the highest possible depth

           `gauss-chebyshev1`  `gauss-chebyshev1-odd`
           `gauss-chebyshev2`  `gauss-chebyshev2-odd`
             approximation using roots of polynomials orthogonal in
             measures  1/sqrt(1-x^2) and sqrt(1-x^2)  (respectively)

          `gauss-gegenbauer`  `gauss-gegenbauer-odd`
            approximation using roots of polynomials orthogonal in
            measure (1-x^2)^alpha

          `gauss-jacobi`
            approximation using roots of polynomials orthogonal in
            measure (1-x)^alpha * (1+x)^beta

          `gauss-laguerre`
            approximation using roots of polynomials orthogonal in
            measure x^alpha * epx(-x)

          `gauss-hermite`  `gauss-hermite-odd`
            approximation using roots of polynomials orthogonal in
            measure |x|^alpha * epx(-x^2)

anisotropic_weights: list or array of weights (Int)
                     length must be dimension or 2*dimension
                     the first dimension weights must be positive
                     see the manual for details

alpha, beta: Float64
      alpha : the alpha parameter for Gegenbauer, Jacobi,
               Hermite and Laguerre rules
      beta  : the beta parameter for Jacobi rules

custom_filename: string giving the path to the file with
             custom-tabulated rule
"""
function makeGlobalGrid(; dimension, outputs, depth, type, rule, anisotropic_weights=[], alpha=0.0, beta=0.0, custom_filename="", level_limits=[])
    tsg = TasmanianSG()
    makeGlobalGrid!(tsg, dimension = dimension, outputs = outputs, depth = depth, type = type, rule = rule , anisotropic_weights = anisotropic_weights,
                    alpha = alpha, beta = beta, custom_filename = custom_filename, level_limits = level_limits)
    return tsg
end

"""
    makeGlobalGrid!(tsg::TasmanianSG; dimension, outputs, depth, type, rule, anisotropic_weights=[], alpha=0.0, beta=0.0, custom_filename="", level_limits=[])

creates a new sparse grid using a global polynomial rule
discards any existing grid held by `tsg`

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

type: string identifying the tensor selection strategy
     `level`     `curved`     `hyperbolic`     `tensor`
     `iptotal`   `ipcurved`   `iphyperbolic`   `iptensor`
     `qptotal`   `qpcurved`   `qphyperbolic`   `qptensor`

type: string identifying the tensor selection strategy
     `level`     `curved`     `hyperbolic`     `tensor`
     `iptotal`   `ipcurved`   `iphyperbolic`   `iptensor`
     `qptotal`   `qpcurved`   `qphyperbolic`   `qptensor`

rule: string (defines the 1-D rule that induces the grid)

       Interpolation rules

          Note: the quadrature induced by those rules is constructed
                by integrating the interpolant

          `clenshaw-curtis`    `clenshaw-curtis-zero`      `fejer2`
          `rleja`    `rleja-odd`  `rleja-double2`   `rleja-double4`
          `rleja-shifted`   `rleja-shifted-even`
          `max-lebesgue`    `max-lebesgue-odd`
          `min-lebesgue`    `min-lebesgue-odd`
          `leja`            `leja-odd`
          `min-delta`       `min-delta-odd`

          `chebyshev`       `chebyshev-odd`
            approximation using roots of Chebyshev polynomials
            non-nested case (in contrast to Clenshaw-Curtis nodes)

       Quadrature rules, the weights target exactness with respect
                         to the highest polynomial degree possible

           `gauss-legendre`  `gauss-legendre-odd`
            approximation using roots of polynomials orthogonal in
            measure Uniform

           `gauss-patterson`  (a.k.a. nested Gauss-Legendre)
            Note: the nodes and weights are hard-coded hence there
            is a limit on the highest possible depth

           `gauss-chebyshev1`  `gauss-chebyshev1-odd`
           `gauss-chebyshev2`  `gauss-chebyshev2-odd`
             approximation using roots of polynomials orthogonal in
             measures  1/sqrt(1-x^2) and sqrt(1-x^2)  (respectively)

          `gauss-gegenbauer`  `gauss-gegenbauer-odd`
            approximation using roots of polynomials orthogonal in
            measure (1-x^2)^alpha

          `gauss-jacobi`
            approximation using roots of polynomials orthogonal in
            measure (1-x)^alpha * (1+x)^beta

          `gauss-laguerre`
            approximation using roots of polynomials orthogonal in
            measure x^alpha * epx(-x)

          `gauss-hermite`  `gauss-hermite-odd`
            approximation using roots of polynomials orthogonal in
            measure |x|^alpha * epx(-x^2)

anisotropic_weights: list or array of weights (Int)
                     length must be dimension or 2*dimension
                     the first dimension weights must be positive
                     see the manual for details

alpha, beta: Float64
      alpha : the alpha parameter for Gegenbauer, Jacobi,
               Hermite and Laguerre rules
      beta  : the beta parameter for Jacobi rules

custom_filename: string giving the path to the file with
             custom-tabulated rule
"""
function makeGlobalGrid!(tsg; dimension, outputs, depth, type, rule, anisotropic_weights=[],
                         alpha=0.0, beta=0.0, custom_filename="", level_limits=[])
    # If dimension, outputs or depth is missing, the value in tsg is taken instead
    if dimension  <= 0
        throw(TasmanianInputError("ERROR: dimension should be a positive integer"))
    end
    if outputs  < 0  
        throw(TasmanianInputError("ERROR: outputs should be a non-negative integer"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
    if !(type in GlobalTypes)
        throw(TasmanianInputError("ERROR: invalid type, see TasmanianSG.GlobalTypes for list of accepted types"))
    end
    if !(rule in GlobalRules)
        throw(TasmanianInputError("ERROR: invalid global rule, see TasmanianSG.lsTsgGlobalRules for list of accepted global rules"))
    end
    panisotropic_weights = check_anisotropic_weights_(anisotropic_weights, dimension, type)
    plevel_limits = check_level_limits_(level_limits, dimension)
    tsgMakeGlobalGrid(tsg.pGrid, dimension, outputs, depth, type, rule, panisotropic_weights, alpha, beta, custom_filename, plevel_limits)
    tsg.dimension = dimension
    tsg.outputs = outputs
    tsg.depth = depth
    return nothing
end

"""
    makeSequenceGrid(dimension, outputs, depth, type, rule, anisotropic_weights=[], level_limits=[])

creates a new sparse grid using a sequence rule
discards any existing grid held by this class

type: string identifying the tensor selection strategy
      'level'     'curved'     'hyperbolic'     'tensor'
      'iptotal'   'ipcurved'   'iphyperbolic'   'iptensor'
      'qptotal'   'qpcurved'   'qphyperbolic'   'qptensor'

rule: string (defines the 1-D rule that induces the grid)
      'leja'       'rleja'      'rleja-shifted'
      'max-lebesgue'   'min-lebesgue'   'min-delta'

anisotropic_weights: list or array of weights
                     length must be `tsg.dimension` or `2*tsg.dimension`
                     the first `tsg.dimension` weights must be positive
                     see the manual for details
"""
function makeSequenceGrid(; dimension, outputs, depth, type, rule, anisotropic_weights=[], level_limits=[])
    tsg = TasmanianSG()
    makeSequenceGrid!(tsg, dimension = dimension, outputs = outputs, depth = depth, type = type,
                      rule = rule, anisotropic_weights = anisotropic_weights, level_limits = level_limits)
    return tsg
end

"""
    makeSequenceGrid!(grid; dimension, outputs, depth, type, rule, anisotropic_weights=[], level_limits=[])

creates a new sparse grid using a sequence rule
discards any existing grid held by this class

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

type: string identifying the tensor selection strategy
      'level'     'curved'     'hyperbolic'     'tensor'
      'iptotal'   'ipcurved'   'iphyperbolic'   'iptensor'
      'qptotal'   'qpcurved'   'qphyperbolic'   'qptensor'

rule: string (defines the 1-D rule that induces the grid)
      'leja'       'rleja'      'rleja-shifted'
      'max-lebesgue'   'min-lebesgue'   'min-delta'

anisotropic_weights: list or array of weights
                     length must be `tsg.dimension` or `2*tsg.dimension`
                     the first `tsg.dimension` weights must be positive
                     see the manual for details
"""
function makeSequenceGrid!(tsg; dimension, outputs, depth, type, rule, anisotropic_weights=[], level_limits=[])
    if dimension <= 0
        throw(TasmanianInputError("ERROR: dimension should be a positive integer"))
    end
    if outputs < 0
        throw(TasmanianInputError("ERROR: outputs should be a non-negative integer"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
    if !(type in GlobalTypes)
        throw(TasmanianInputError("ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types"))
    end
    if !(rule in SequenceRules)
        throw(TasmanianInputError("ERROR: invalid sequence rule, see TasmanianSG.lsTsgSequenceRules for list of accepted sequence rules"))
    end

    panisotropic_weights = check_anisotropic_weights_(anisotropic_weights, dimension, type)
    plevel_limits = check_level_limits_(level_limits, dimension)

    tsgMakeSequenceGrid(tsg.pGrid, dimension, outputs, depth, type, rule, panisotropic_weights, plevel_limits)
    tsg.dimension = dimension
    tsg.outputs = outputs
    tsg.depth = depth
    return nothing
end

"""
    makeLocalPolynomialGrid(; dimension, outputs, depth, order=1, rule="localp", level_limits=[])

creates a new sparse grid using a local polynomial rule
discards any existing grid held by TSG

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

order: int (must be -1 or bigger)
        -1 indicates largest possible order
         1 means linear, 2 means quadratic, etc.
         0 means piece-wise constant, it has different hierarchy
           then the other orders, most notably the 1D rule
           triples the number of points per level (as opposed
           to double for the other cases)

rule: string (defines the 1-D rule that induces the grid)
      `localp` `localp-zero`  `semi-localp`  `localp-boundary`

"""
function makeLocalPolynomialGrid(; dimension, outputs, depth, order=1, rule="localp", level_limits=[])
    tsg = TasmanianSG()
    makeLocalPolynomialGrid!(tsg, dimension = dimension, outputs = outputs, depth = depth, order=order, rule=rule, level_limits=level_limits)
    return tsg
end

"""
    makeLocalPolynomialGrid!(tsg::TasmanianSG; dimension, outputs, depth, order=1, rule="localp", level_limits=[])

creates a new sparse grid using a local polynomial rule
discards any existing grid held by TSG

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

order: int (must be -1 or bigger)
        -1 indicates largest possible order
         1 means linear, 2 means quadratic, etc.
         0 means piece-wise constant, it has different hierarchy
           then the other orders, most notably the 1D rule
           triples the number of points per level (as opposed
           to double for the other cases)

rule: string (defines the 1-D rule that induces the grid)
      `localp` `localp-zero`  `semi-localp`  `localp-boundary`

"""
function makeLocalPolynomialGrid!(tsg::TasmanianSG; dimension, outputs, depth, order=1, rule="localp", level_limits=[])
    if dimension <= 0
        throw(TasmanianInputError("ERROR: dimension should be a positive integer"))
    end
    if outputs < 0
        throw(TasmanianInputError("ERROR: outputs should be a non-negative integer"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
    if order < -1
        throw(TasmanianInputError("ERROR: order should be a non-negative integer"))
    end
    if !(rule in LocalRules)
        throw(TasmanianInputError("ERROR: invalid local polynomial rule, see TasmanianSG.LocalRules for list of accepted sequence rules"))
    end
    plevel_limits = check_level_limits_(level_limits, dimension)
    tsgMakeLocalPolynomialGrid(tsg.pGrid, dimension, outputs, depth, order, rule, plevel_limits)
    tsg.dimension = dimension
    tsg.outputs = outputs
    tsg.depth = depth
    return nothing
end

"""
    makeWaveletGrid(; dimension, outputs, depth, order=1, level_limits=[])

creates a new sparse grid using a wavelet rule
discards any existing grid held by `tsg`

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

order: Int (must be 1 or 3)
       only wavelets of order 1 and 3 are implemented
"""
function makeWaveletGrid(; dimension, outputs, depth, order=1, level_limits=[])
    tsg = TasmanianSG()
    makeWaveletGrid!(tsg, dimension = dimension, outputs = outputs, depth = depth, order=order, level_limits=level_limits)
    return tsg
end

"""
    makeWaveletGrid!(tsg::TasmanianSG; dimension, outputs, depth, order=1, level_limits=[])

creates a new sparse grid using a wavelet rule
discards any existing grid held by `tsg`

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

order: Int (must be 1 or 3)
       only wavelets of order 1 and 3 are implemented
"""
function makeWaveletGrid!(tsg; dimension, outputs, depth, order=1, level_limits=[])
    if dimension <= 0
        throw(TasmanianInputError("ERROR: dimension should be a positive integer"))
    end
    if outputs < 0
        throw(TasmanianInputError("ERROR: outputs should be a non-negative integer"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
    if !(order in [1, 3])
        throw(TasmanianInputError("ERROR: order should be either 1 or 3 (only linear and cubic wavelets are available)"))
    end
    
    plevel_limits = check_level_limits_(level_limits, dimension)

    tsgMakeWaveletGrid(tsg.pGrid, dimension, outputs, depth, order, plevel_limits)
    tsg.dimension = dimension
    tsg.outputs = outputs
    tsg.depth = depth
    return nothing
end

"""
    makeFourierGrid(; dimension, outputs, depth, type, anisotropic_weights=[], level_limits=[])

creates a new sparse grid using a Fourier rule
discards any existing grid held by this class

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

type: string identifying the tensor selection strategy
      'level'     'curved'     'hyperbolic'     'tensor'
      'iptotal'   'ipcurved'   'iphyperbolic'   'iptensor'
      'qptotal'   'qpcurved'   'qphyperbolic'   'qptensor'

anisotropic_weights: list or array of weights
                     length must be iDimension or 2*iDimension
                     the first tsg.dimension weights must be positive
                     see the manual for details
"""
function makeFourierGrid(; dimension, outputs, depth, type, anisotropic_weights=[], level_limits=[])
    tsg = TasmanianSG()
    makeFourierGrid!(tsg, dimension = dimension, outputs = outputs, depth = depth, type = type,
                     anisotropic_weights = anisotropic_weights, level_limits = level_limits)
    return tsg
end

"""
    makeFourierGrid!(tsg::TasmanianSG; dimension, outputs, depth, type, anisotropic_weights=[], level_limits=[])

creates a new sparse grid using a Fourier rule
discards any existing grid held by this class

dimension: int (positive)
           the number of inputs

outputs: int (non-negative)
         the number of outputs

depth: int (non-negative)
       controls the density of the grid, i.e.,
       the offset for the tensor selection, the meaning of
       depth depends on type
       Example 1: type == 'iptotal' will give a grid that
                  interpolates exactly all polynomials of
                  degree up to and including depth
       Example 2: type == 'qptotal' will give a grid that
                  integrates exactly all polynomials of degree
                  up to and including depth

type: string identifying the tensor selection strategy
      'level'     'curved'     'hyperbolic'     'tensor'
      'iptotal'   'ipcurved'   'iphyperbolic'   'iptensor'
      'qptotal'   'qpcurved'   'qphyperbolic'   'qptensor'

anisotropic_weights: list or array of weights
                     length must be iDimension or 2*iDimension
                     the first tsg.dimension weights must be positive
                     see the manual for details
"""
function makeFourierGrid!(tsg::TasmanianSG; dimension, outputs, depth, type, anisotropic_weights=[], level_limits=[])
    if dimension <= 0
        throw(TasmanianInputError("ERROR: dimension should be a positive integer"))
    end
    if outputs < 0
        throw(TasmanianInputError("ERROR: outputs should be a non-negative integer"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
    if !(type in GlobalTypes)
        throw(TasmanianInputError("ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types"))
    end
    panisotropic_weights = check_anisotropic_weights_(anisotropic_weights, dimension, type)
    plevel_limits = check_level_limits_(level_limits, dimension)
    tsgMakeFourierGrid(tsg.pGrid, dimension, outputs, depth, type, panisotropic_weights, plevel_limits)
    tsg.dimension = dimension
    tsg.outputs = outputs
    tsg.depth = depth
    return tsg
end

"""
    copyGrid(tsg::TasmanianSG, outputs_begin = 0, outputs_end = -1)

accepts an instance of TasmanianSparseGrid class and creates
a hard copy of the grid and all included data
original grid is not modified

tsg: instance of TasmanianSG
     the source for the copy

outputs_begin: integer indicating the first output to copy

outputs_end: integer one bigger than the last output to copy
            if set to -1, all outputs from outputs_begin to
            the end will be copied

Examples:

newgrid = copyGrid(other, 0, -1) # copy all outputs (default)
newgrid = copyGrid(other, 0, getNumOutputs(other)) # also copy all
newgrid = copyGrid(other, 0, 3) # copy outputs 0, 1, and 2
newgrid = copyGrid(other, 1, 4) # copy outputs 1, 2, and 3
"""
function copyGrid(tsg::TasmanianSG, outputs_begin = 0, outputs_end = -1)
    newgrid = TasmanianSG()
    tsgCopySubGrid(newgrid.pGrid, tsg.pGrid, outputs_begin, outputs_end)
    return newgrid
end

"""
    copyGrid!(newgrid::TasmanianSG, oldgrid::TasmanianSG, outputs_begin = 0, outputs_end = -1)

makes a hard copy of old grid on new grid. The old grid is not modified

new: instance of TasmanianSG
     the destination for the copy

tsg: instance of TasmanianSG
     the source for the copy

outputs_begin: integer indicating the first output to copy

outputs_end: integer one bigger than the last output to copy
            if set to -1, all outputs from outputs_begin to
            the end will be copied

Examples:

copyGrid!(newgrid, other, 0, -1) # copy all outputs (default)
copyGrid!(newgrid, other, 0, getNumOutputs(other)) # also copy all
copyGrid!(newgrid, other, 0, 3) # copy outputs 0, 1, and 2
copyGrid!(newgrid, other, 1, 4) # copy outputs 1, 2, and 3
"""
function copyGrid!(new::TasmanianSG, tsg::TasmanianSG, outputs_begin = 0, outputs_end = -1)
    newgrid = TasmanianSG()
    tsgCopySubGrid(newgrid.pGrid, tsg.pGrid, outputs_begin, outputs_end)
    return nothing
end

"""
  updateGlobalGrid!(tsg::TasmanianSG; depth, type, anisotropic_weights=[], level_limits=[])
adds the points defined by depth, type and anisotropy
to the existing grid

basically, the same as calling makeGlobalGrid with rule,
           alpha and beta of this grid and the new depth,
           type and anisotropic_weights then adding the
           resulting points to the current grid

inputs: see help(makeGlobalGrid)
"""
function updateGlobalGrid!(tsg::TasmanianSG; depth, type, anisotropic_weights=[], level_limits=[])
    if !isGlobal(tsg)
        throw(TasmanianInputError("ERROR: calling updateGlobalGrid for a grid that is not global"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
        if !(type in GlobalTypes)
            throw(TasmanianInputError("ERROR: invalid type, see GlobalTypes for list of accepted types"))
        end

    panisotropic_weights = check_anisotropic_weights_(anisotropic_weights,  tsg.dimension, type)
    plevel_limits = check_level_limits_(level_limits, tsg.dimension)

    tsgUpdateGlobalGrid(tsg.pGrid, depth, type, panisotropic_weights, plevel_limits)
end

"""
    updateSequenceGrid!(tsg::TasmanianSG; depth, type, anisotropic_weights=[], level_limits=[])

adds the points defined by depth, type and anisotropy
to the existing grid

basically, the same as calling makeSequenceGrid() with rule,
           of this grid and the new depth, type and
           anisotropic_weights then adding the resulting points
           to the current grid

inputs: see help(makeSequenceGrid)
"""
function updateSequenceGrid!(tsg::TasmanianSG; depth, type, anisotropic_weights=[], level_limits=[])
    if !isSequence(tsg)
        throw(TasmanianInputError("ERROR: calling updateSequenceGrid for a grid that is not a sequence grid"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
    if !(type in GlobalTypes)
        throw(TasmanianInputError("ERROR: invalid type, see GlobalTypes for list of accepted types"))
    end

    panisotropic_weights = check_anisotropic_weights_(anisotropic_weights, tsg.dimension, type)
    plevel_limits = check_level_limits_(level_limits, tsg.dimension)

    tsgUpdateSequenceGrid(tsg.pGrid, depth, type, panisotropic_weights, plevel_limits)
end

"""
    updateFourierGrid!(tsg::TasmanianSG; depth, type, anisotropic_weights=[], level_limits=[])
adds the points defined by depth, type and anisotropy
to the existing grid

basically, the same as calling makeFourierGrid() with rule,
           of this grid and the new depth, type and
           anisotropic_weights then adding the resulting points
           to the current grid

inputs: see help(makeGlobalGrid)
"""
function updateFourierGrid!(tsg::TasmanianSG; depth, type, anisotropic_weights=[], level_limits=[])
    if !isFourier(tsg)
        throw(TasmanianInputError("ERROR: calling updateFourierGrid for a grid that is not a Fourier grid"))
    end
    if depth < 0
        throw(TasmanianInputError("ERROR: depth should be a non-negative integer"))
    end
    if !(type in GlobalTypes)
            throw(TasmanianInputError("ERROR: invalid type, see Tasmanian.TasmanianSG.lsTsgGlobalTypes for list of accepted types"))
    end

    panisotropic_weights = check_anisotropic_weights_(anisotropic_weights, tsg.dimension, type)
    plevel_limits = check_level_limits_(level_limits, tsg.dimension)

    tsgUpdateFourierGrid(tsg.pGrid, depth, type, panisotropic_weights, plevel_limits)
end

"""
    getAlpha(tsg::TasmanianSG)

    returns the value of alpha in the call to makeGlobalGrid
    if makeGlobalGrid has not been called, returns 0.0
"""
function getAlpha(tsg::TasmanianSG)
    return tsgGetAlpha(tsg.pGrid)
end

"""
    getBeta(tsg::TasmanianSG)

    returns the value of beta in the call to makeGlobalGrid
    if makeGlobalGrid has not been called, returns 0.0
"""
function getBeta(tsg::TasmanianSG)
    return tsgGetBeta(tsg.pGrid)
end

"""
    getOrder(tsg::TasmanianSG)

returns the value of iOrder in the call to
makeLocalPolynomialGrid or makeWaveletGrid
if makeLocalPolynomialGrid and makeWaveletGrid
have not been called, returns -1
"""
getOrder(tsg::TasmanianSG) = convert(Int, tsgGetOrder(tsg.pGrid))

"""
    getNumDimensions(tsg::TasmanianSG)

returns the value of iDimension in the make***Grid command
if no grid has been made, it returns 0
"""
getNumDimensions(tsg::TasmanianSG) = convert(Int, tsgGetNumDimensions(tsg.pGrid))

"""
    getNumOutputs(tsg::TasmanianSG)

returns the value of iOutputs in the make***Grid command
if no grid has been made, it returns 0
"""
getNumOutputs(tsg::TasmanianSG) = convert(Int, tsgGetNumOutputs(tsg.pGrid))

"""
    getRule(tsg::TasmanianSG)

returns the value of rule in the make***Grid command
if makeWaveletGrid is used, returns "wavelet"
if no grid has been made, it returns "unknown"
"""
function getRule(tsg::TasmanianSG)
    buffer_size = 128
    name = Vector{Cchar}(undef, buffer_size)
    num_actual = Ref{Cint}()
    tsgCopyRuleChars(tsg.pGrid, buffer_size, name, num_actual)
    return String(Char.(name[1:num_actual[]]))
end

"""
    getNumLoaded(tsg::TasmanianSG)

returns the number of points loaded in the existing interpolant
"""
getNumLoaded(tsg::TasmanianSG) = convert(Int, tsgGetNumLoaded(tsg.pGrid))

"""
    getNumNeeded(tsg::TasmanianSG)

returns the number of points needed to form the interpolant or
        form the next interpolant following a refinement
"""
getNumNeeded(tsg::TasmanianSG) = convert(Int, tsgGetNumNeeded(tsg.pGrid))

"""
    getNumPoints(tsg::TasmanianSG)

if points have been loaded, returns the same as getNumLoaded()
otherwise, returns the same as getNumNeeded()
"""
getNumPoints(tsg::TasmanianSG) = convert(Int,tsgGetNumPoints(tsg.pGrid))

"""
    getLoadedPoints(tsg::TasmanianSG)

returns the points loaded in the existing interpolant

output: matrix of size getNumDimensions() X getNumNeeded()
    each column  corresponds to one point
"""
function getLoadedPoints(tsg::TasmanianSG)
    NumDimensions = getNumDimensions(tsg)
    NumPoints = getNumLoaded(tsg)
    outputs = zeros(NumDimensions, NumPoints)
    NumPoints > 0 && tsgGetLoadedPointsStatic(tsg.pGrid, outputs)
    return outputs
end

"""
    getNeededPoints(tsg::TasmanianSG)

returns the points needed to form the interpolant or the next
level of refinement following a set***Refinement() call

output: 2-D array of size dimension X getNumNeeded()
    each column  corresponds to one point
    if (getNumNeeded() == 0): returns zeros(dimensions, 0)
"""
function getNeededPoints(tsg::TasmanianSG)
    NumDimensions = getNumDimensions(tsg)
    NumPoints = getNumNeeded(tsg)
    outputs = zeros(NumDimensions, NumPoints)
    NumPoints > 0 && tsgGetNeededPointsStatic(tsg.pGrid, outputs)
    return outputs
end

"""
    getPoints(tsg::TasmanianSG)

if points have been loaded, gives the same as getLoadedPoints()
otherwise, returns the same as getNeededPoints()
"""
function getPoints(tsg::TasmanianSG)
    NumDimensions = getNumDimensions(tsg)
    NumPoints = getNumPoints(tsg)
    outputs = zeros(NumDimensions, NumPoints)
    NumPoints > 0 && tsgGetPointsStatic(tsg.pGrid, outputs)
    return outputs
end

"""
    getQuadratureWeights(tsg)

returns the quadrature weights associated with
the points in getPoints()

output: a vector of length getNumPoints()
        the order of the weights matches
        the order in getPoints()
"""
function getQuadratureWeights(tsg)
    NumPoints = getNumPoints(tsg)
    weights = zeros(NumPoints)
    NumPoints  > 0 && tsgGetQuadratureWeightsStatic(tsg.pGrid, weights)
    return weights
end
    
"""
    getInterpolationWeights(tsg, x)

returns the interpolation weights associated with the points
in getPoints()

x: a vector with length dimensions
   the entries indicate the points for evaluating the weights

output: a vector of length getNumPoints()
        the order of the weights matches the order in getPoints()
"""
function getInterpolationWeights(tsg, x)
    NumX = length(x)
    if NumX != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: length(x) should equal $(getNumDimensions(tsg)) instead it equals $NumX"))
    end
    NumPoints = getNumPoints(tsg)
    weights = zeros(NumPoints)
    NumPoints > 0 && tsgGetInterpolationWeightsStatic(tsg.pGrid, x, weights)
    return weights
end

"""
    getInterpolationWeightsBatch(tsg, x):
returns the interpolation weights associated with the points
in getPoints()

finds multiple weights with a single library call
uses OpenMP if enabled in libtasmaniansparsegrids.so

x: a matrix with first dimension getNumDimensions()
   each column in the array is a single requested point

output: a matrix
        with dimensions getNumPoints() X size(x, 2)  
        each row corresponds to the weight for one row of x
"""
function getInterpolationWeightsBatch(tsg, x)
    if ndims(x) > 2
        throw(TasmanianInputError("ERROR: x should be a vector or a matrix instead it has dimension $(ndims(x))"))
    end
    NumDim = size(x, 1)
    NumX = size(x, 2)
    if NumDim != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: size(x, 1) should equal $(getNumDimensions(tsg)) instead it equals $NumDim"))
    end
    weights = zeros(getNumPoints(tsg), NumX)
    NumDim > 0 && NumX > 0 && tsgBatchGetInterpolationWeightsStatic(tsg.pGrid, x, NumX, weights)
    return weights
end

"""
    loadNeededPoints!(tsg::TasmanianSG, vals::Array{Float64})

loads the values of the target function at the needed points
if there are no needed points, this reset the currently loaded
values

vals: an array with dimensions outputs X getNumNeeded() 
      each column corresponds to the values of the outputs at
      the corresponding needed point. The order and leading
      dimension must match the points obtained from
      getNeededPoints()
"""
function loadNeededPoints!(tsg::TasmanianSG, vals::AbstractArray{Float64})
    numOutputs = getNumOutputs(tsg)
    nd = ndims(vals)
    if nd == 1
        n = length(vals)
        if mod(n, numOutputs) == 0
            n1 = numOutputs
            n2 = n/numOutputs
        else
            throw(TasmanianInputError("ERROR: vals is a vector but its length isn't a multiple of the numOutputs"))
        end
    elseif nd == 2
        n1, n2 = size(vals)
    else
        throw(TasmanianInputError("ERROR: vals must be a vector or a matrix"))
    end
    if n1 != numOutputs
        throw(TasmanianInputError("ERROR: leading dimension of vals is $n1 but the number of outputs is set to $(getNumOutputs(tsg))"))
    end
    if getNumNeeded(tsg) == 0
        if n2 != getNumLoaded(tsg)
            throw(TasmanianInputError("ERROR: the second dimension of vals is $n2 but the number of current points is $(getNumLoaded(tsg))"))
        end
    elseif n2 != getNumNeeded(tsg)
        throw(TasmanianInputError("ERROR: the second dimension of vals is $n2 but the number of needed points is $(getNumNeeded(tsg))"))
    end
    tsgLoadNeededPoints(tsg.pGrid, vals)
end

"""
    loadNeededValues!(tsg::TasmanianSG, vals::Array{Float64})

Alias of loadNeededPoints()
"""
loadNeededValues!(tsg::TasmanianSG, vals::Array{Float64}) = loadNeededPoints!(tsg, vals)

"""
    getLoadedValues(tsg::TasmanianSG)
Returns the model values as given to Tasmanian by the loadNeededPoints() method.
The ordering will match the current internal ordering, e.g., mixing the different
model values from different refinement iterations.

Returns a matrix with size getNumOutputs() by getNumPoints()
"""
function getLoadedValues(tsg::TasmanianSG)
    NumPoints  = getNumPoints(tsg)
    NumOutputs = getNumOutputs(tsg)
    vals = zeros(NumOutputs, NumPoints)
    NumPoints > 0 && NumOutputs > 0 && tsgGetLoadedValuesStatic(tsg.pGrid, vals)
    return vals
end

"""
evaluateThreadSafe(tsg::TasmanianSG, x)
        evaluates the intepolant at a single points of interest and
        returns the result
        This is the thread safe version, but it does not use
        acceleration of any type

        this should be called after the grid has been created and after
        values have been loaded

        x: vector with length dimensions
           the entries indicate the points for evaluating the weights

        output: returns a vector of length getNumOutputs
                the values of the interpolant at x
"""
function evaluateThreadSafe(tsg::TasmanianSG, x)
    if (getNumLoaded(tsg) == 0)
        throw(TasmanianInputError("ERROR: cannot call evaluate for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    if ndims(x) != 1
        throw(TasmanianInputError("ERROR: x should be a vector"))
    end
    NumX = length(x)
    if NumX != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: x should have length $(getNumDimensions(tsg)) instead it has length $NumX"))
    end
    NumOutputs = getNumOutputs(tsg)
    y = zeros(NumOutputs)
    NumOutputs > 0 && tsgEvaluate(tsg.pGrid, x, y)
    return y 
end

"""
    evaluate(tsg::TasmanianSG, x)
evaluates the intepolant at a single points of interest and
returns the result
This is the accelerated version using the selected acceleration
type, but it is potentially not thread safe

this should be called after the grid has been created and after
values have been loaded

x: a vector with length getNumDimensions()
   the entries indicate the points for evaluating the weights

output: returns vector of length getNumOutputs()
        the values of the interpolant at fX
"""
function evaluate(tsg::TasmanianSG, x)
    if getNumLoaded(tsg) == 0
        throw(TasmanianInputError("ERROR: cannot call evaluate for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    if ndims(x) != 1
        throw(TasmanianInputError("ERROR: x should be a vector"))
    end
    NumX = length(x)
    if NumX != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: x should have length $(getNumDimensions(tsg)) instead it has length $NumX"))
    end
    NumOutputs = getNumOutputs(tsg)
    y = Vector{Float64}(undef, NumOutputs)
    NumOutputs > 0 && NumX > 0 && tsgEvaluateFast(tsg.pGrid, x, y)
    return y
end

"""
     evaluateBatch(tsg::TasmanianSG, vals)

evaluates the intepolant at the points of interest and returns
the result

this should be called after the grid has been created and after
values have been loaded

vals: a vector or a matrix
      with first dimension equal to dimensions
      each column in the array is a single requested point

output: a vector or a matrix
        with dimensions outputs X size(vals, 2) 
        each columns corresponds to the value of the interpolant
        for one columns of vals
"""
function evaluateBatch(tsg::TasmanianSG, vals)
    NumOutputs = getNumOutputs(tsg)
    NumX = size(vals, 2)
    y = Matrix{Float64}(undef, NumOutputs, NumX)
    evaluateBatch!(y, tsg, vals)
    return y
end

"""
     evaluateBatch!(y::AbstractVecOrMat{Float64}, tsg::TasmanianSG, vals::AbstractVecOrMat{Float64})

evaluates the intepolant at the points of interest and set the result in `y`

this should be called after the grid has been created and after
values have been loaded

y: a vector or a matrix
   with dimensions outputs X size(vals, 2) 
   each columns corresponds to the value of the interpolant
   for one columns of vals

vals: a vector or a matrix
      with first dimension equal to dimensions
      each column in the array is a single requested point
"""
function evaluateBatch!(y::AbstractVecOrMat{Float64}, tsg::TasmanianSG, vals::AbstractVecOrMat{Float64})
    !isa(y, StridedArray) && throw(TasmanianInputError("ERROR: y must be a StridedArray"))
    !isa(vals, StridedArray) && throw(TasmanianInputError("ERROR: vals must be a StridedArray"))
    if !isa(vals, AbstractVecOrMat{Float64})
#        vals = convert(AbstractVecOrMat{Float64}, vals)
    end
    if getNumLoaded(tsg) == 0
        throw(TasmanianInputError("ERROR: cannot call evaluateBatch for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    n1 = size(vals)
    if ndims(vals) > 2
        throw(TasmanianInputError("ERROR: vals should be a vector or a matrix instead it has $(ndims(vals)) dimensions"))
    end
    if ndims(y) != ndims(vals)
        throw(TasmanianInputError("ERROR: y and vals must have the same number of dimensions"))
    end
    NumDim = size(vals, 1)
    NumX = (ndims(vals) == 2) ? size(vals, 2) : 1
    NumDimOut = size(y, 1)
    NumXOut = (ndims(y) == 2) ? size(y, 2) : 1
    if NumX > NumXOut
        throw(TasmanianInputError("ERROR: y must have at least as many columns as vals"))
    end
    if (NumDim != getNumDimensions(tsg))
        throw(TasmanianInputError("ERROR: size(vals, 1) should equal $(getNumDimensions(tsg)) instead it equals $NumDim"))
    end
    if (NumDimOut != getNumOutputs(tsg))
        throw(TasmanianInputError("ERROR: size(y, 1) should equal $(getNumOutputs(tsg)) instead it equals $NumDimOut"))
    end
    NumX > 0 && NumDim > 0 && tsgEvaluateBatch( tsg.pGrid, vals, NumX, y)
    return y
end

"""
    integrate(tsg::TasmanianSG)

returns the integral of the interpolant

output: returns a vector of length outputs
        the integral of the interpolant
"""
function integrate(tsg::TasmanianSG)
    if getNumLoaded(tsg) == 0
        throw(TasmanianInputError("ERROR: cannot call integrate for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    NumOutputs = getNumOutputs(tsg)
    aQ = zeros(NumOutputs)
    tsgIntegrate(tsg.pGrid, aQ)
    return aQ
end

"""
    differentiate!(jacobian::VecOrMat{Float64}, tsg::TasmanianSG, x::AbstractVector{Float64})

returns the derivative (Jacobian or gradient vector) of the interpolant

jacobian:  a vector or a matrix
           with dimensions dimensions X outputs
           each row corresponds to the value of the interpolant
           for one columns of vals
tsg:  an instance of TasmanianSG
x: a vector with length dimensions which is the evaluation point
"""
function differentiate!(jacobian::VecOrMat{Float64}, tsg::TasmanianSG, x::Vector{Float64})
    dimensions = getNumDimensions(tsg)
    outputs = getNumOutputs(tsg)
    if dimensions != size(jacobian, 1)
        throw(TasmanianInputError("ERROR: jacobian must have as many rows as getNumDimensions(tsg)"))
    end
    if outputs != size(jacobian, 2)
        throw(TasmanianInputError("ERROR: jacobian must have as many columns as getNumOutputs(tsg)"))
    end
    if dimensions != length(x)
        throw(TasmanianInputError("ERROR: x must have length equal to getNumDimensions(tsg)"))
    end
    tsgDifferentiate( tsg.pGrid, x, jacobian)
    return jacobian
end

"""
    differentiate(tsg::TasmanianSG, x::Vector{Float64})

returns the derivative (Jacobian or gradient vector) of the interpolant
as a vector or a matrix with dimensions dimensions X outputs
each row corresponds to the value of the interpolant

tsg:  an instance of TasmanianSG
x: a vector with length iDimensions which is the evaluation point
"""
function differentiate(tsg::TasmanianSG, x::Vector{Float64})
    dimensions = getNumDimensions(tsg)
    outputs = getNumOutputs(tsg)
    if dimensions != length(x)
        throw(TasmaninaInputError("ERROR: x must have as many columns as getNumDimensions(tsg)"))
    end
    jacobian = zeros(dimensions, outputs)
    differentiate!(jacobian, tsg, x)
    return jacobian
end

"""
    isGlobal(tsg::TasmanianSG)

returns true if using a global grid
"""
isGlobal(tsg::TasmanianSG)          = convert(Bool, tsgIsGlobal(tsg.pGrid))

"""
    isSequence(tsg::TasmanianSG)

returns true if using a sequence grid
"""
isSequence(tsg::TasmanianSG)        = convert(Bool, tsgIsSequence(tsg.pGrid))

"""
    isLocalPolynomial(tsg::TasmanianSG)

returns true if using a local polynomial grid
"""
isLocalPolynomial(tsg::TasmanianSG) = convert(Bool,tsgIsLocalPolynomial(tsg.pGrid))

"""
    isWavelet(tsg::TasmanianSG)
returns true if using a local wavelet grid
"""
isWavelet(tsg::TasmanianSG)         = convert(Bool,tsgIsWavelet(tsg.pGrid))

"""
    isFourier(tsg::TasmanianSG)

returns true if using a Fourier grid
"""
isFourier(tsg::TasmanianSG)         = convert(Bool,tsgIsFourier(tsg.pGrid))

"""
    setDomainTransform!(tsg::TasmanianSG, Transform::VecOrMat)

sets the lower and upper bound for each dimension

Note: gauss-laguerre and gauss-hermite rules are defined on
      unbounded domain, in which case  this sets the
      shift and scale parameters, consult the manual

Transform: a matrix of size dimension X 2
           transform specifies the lower and upper bound
           of the domain in each direction.

           For gauss-laguerre and gauss-hermite grids, the
           transform gives the a and b parameters of the
           weights
            exp(-b (x - a))
            exp(-b (x - a)^2)
"""
function setDomainTransform!(tsg::TasmanianSG, transformation::VecOrMat)
    n = size(transformation)
    if length(n) != 2
        throw(TasmanianInputError("ERROR: transformation should be a matrix"))
    end
    if n[1] != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: the first dimension of transformation is $(n[1]) and it should match iDimension: $(getNumDimensions(tsg))"))
    end
    if n[2] != 2
        throw(TasmanianInputError("ERROR: the second dimension of transformation is $(n[2]) and it should be 2"))
    end
    a = transformation
    b = view(transformation, :, 2)
    tsgSetDomainTransform(tsg.pGrid, a, b)
end

"""
isSetDomainTransform(tsg::TasmanianSG)
returns true if the grid is defined for non-canonical domain
        returns false if using a canonical domain

"""
isSetDomainTransform(tsg::TasmanianSG) = convert(Bool, tsgIsSetDomainTransfrom(tsg.pGrid))

"""
    clearDomainTransform!(tsg::TasmanianSG)
resets the domain to canonical
loaded values will be kept, however, the values now correspond
to canonical points and may be invalid for your application
"""
clearDomainTransform!(tsg::TasmanianSG) = tsgClearDomainTransform(tsg.pGrid)

"""
    getDomainTransform(tsg::TasmanianSG)
returns Transform from the call to setDomainTransform()

if setDomainTransform() has not been called or if the
transformed has been cleared by clearDomainTransform(),
then this returns an empty matrix
"""
function getDomainTransform(tsg::TasmanianSG)
    if !isSetDomainTransform(tsg)
        return zeros(0,2)
    end
    NumDimensions = getNumDimensions(tsg)
    transformation = zeros(NumDimensions, 2)
    a = transformation
    b = view(transformation, :, 2)
    tsgGetDomainTransform(tsg.pGrid, a, b)
    return transformation
end

"""
    setConformalTransformASIN!(tsg::TasmanianSG, truncation)
sets conformal domain transform based on truncated
Maclaurin series of arcsin()

truncation: a vector of non-negative integers
            indicating the truncation order in each direction
            0 indicates no transform applied to this direction
"""
function setConformalTransformASIN!(tsg::TasmanianSG, truncation)
    if ndims(truncation) != 1
        throw(TasmanianInputError("ERROR: truncation should be a vector"))
    end
    if length(truncation) != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: the length of truncation is $(length(truncation)) and it should match dimension: $(getNumDimensions(tsg))"))
    end
    if !isa(truncation, AbstractVector{Float64})
#        truncation = convert(Vector{Float64}, truncation)
    end
    tsgSetConformalTransformASIN(tsg.pGrid, truncation)
end

"""
    isSetConformalTransformASIN(tsg::TasmanianSG)

returns true if conformal transform is set
returns false otherwise

see: setConformalTransformASIN()
"""
isSetConformalTransformASIN(tsg::TasmanianSG) = convert(Bool, tsgIsSetConformalTransformASIN(tsg.pGrid))

"""
    clearConformalTransform!(tsg::TasmanianSG)
resets the conformal domain transform
loaded values will be kept, however, the values now correspond
to canonical points and may be invalid for your application
"""
clearConformalTransform!(tsg::TasmanianSG) = tsgClearConformalTransform(tsg.pGrid)

"""
    getConformalTransformASIN(tsg::TasmanianSG)
returns Truncation from the call to setConformalTransformASIN()

if setConformalTransformASIN() has not been called or if the
transformed has been cleared by clearConformalTransform(),
then this returns an empty vector
"""
function getConformalTransformASIN(tsg::TasmanianSG)
        if !isSetConformalTransformASIN(tsg)
            return Vector{Int32}(undef, 0)
        end
        NumDimensions = getNumDimensions(tsg)
        truncation = Vector{Int32}(undef, NumDimensions)
        tsgGetConformalTransformASIN(tsg.pGrid, truncation)
        return truncation
end

"""
    clearLevelLimits!(tsg::TasmanianSG)

clears the limits set by the last make***Grid or refine command
if no limits are set, this has no effect
"""
clearLevelLimits!(tsg::TasmanianSG) = tsgClearLevelLimits(tsg.pGrid)

"""
    getLevelLimits(tsg::TasmanianSG)

returns the limits set by the last call to make***Grid or refine
returns a vector of integers corresponding to the limits for each
direction, -1 indicates no limit
"""
function getLevelLimits(tsg::TasmanianSG)
    NumDimensions = getNumDimensions(tsg)
    truncation = Vector{Int32}(undef, NumDimensions)
    tsgGetLevelLimits(tsg.pGrid, truncation)
    return truncation
end

"""
    setAnbisotropicRefinement!(tsg::TasmanianSG; type, min_growth, output, level_limits = [])

estimates anisotropic coefficients from the current set of
loaded points and updates the grid with the best points
according to the estimate

type: string identifying the estimate to use (see the Manual)
      recommended: 'iptotal'   'ipcurved'

min_growth: int (positive)
            minimum number of new points to include in the new grid

output: int (indicates the output to use) selects which output to use for refinement
        sequence grids, using -1 indicates to use all outputs

level_limits: (if not empty) will be used to overwrite the currently set limits. The limits must be either empty
              or have size getNumDimensions(); if empty, the current set of limits will be used.
"""
function setAnisotropicRefinement!(tsg::TasmanianSG; type, min_growth, output, level_limits = [])
    if getNumOutputs(tsg) == 0
             throw(TasmanianInputError("ERROR: cannot set refinement for grid with outputs = 0"))
    end
    if getNumLoaded(tsg) == 0
        throw(TasmanianInputError("ERROR: cannot call setAnisotropicRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    if min_growth <= 0
        throw(TasmanianInputError("ERROR: the number min_growth should be positive integer"))
    end
    if output == -1
        if !isSequence(tsg)
            throw(TasmanianInputError("ERROR: output = -1 can be used only for sequence grids"))
        end
    end
    if output < -1
        throw(TasmanianInputError("ERROR: output should be -1 or a non-negative integer"))
    end
        
    if output >= getNumOutputs(tsg)
        throw(TasmanianInputError("ERROR: output cannot exceed the index of the last output $(getNumOutputs(tsg))"))
    end
    if !(type in GlobalTypes)
        throw(TasmanianInputError("ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types"))
    end

    plevel_limits = check_level_limits_(level_limits, tsg.dimension)

    tsgSetAnisotropicRefinement(tsg.pGrid, type, min_growth, output, plevel_limits)
end
             
"""
    getAnisotropicRefinement!(tsg::TasmanianSG; type, min_growth, output, level_Limits = [])

Calls setAnistropicRefinement() on the inputs and then getNeededPoints().
"""
function getAnisotropicRefinement!(tsg::TasmanianSG; type, min_growth, output, level_limits = [])
    setAnisotropicRefinement!(tsg, type, min_growth, output, level_limits=level_limits)
    return getNeededPoints(tsg)
end

"""
    estimateAnisotropicCoefficients(tsg::TasmanianSG; type, output)
returns the estimate of the anisotropic coefficients from the
current set of loaded points
see the manual

type: string identifying the estimate to use (see the Manual)
       recommended: 'iptotal'   'ipcurved'


output: int (indicates the output to use)
     selects which output to use for refinement
     sequence grids accept -1 to indicate all outputs

returns vector of length getNumDimensions() or 2*getNumDimensions()
        the first set of getNumDimensions() entries correspond to the xi coefficients
         the second set of getNumDimensions() entries correspond to the eta coefficients
"""
function estimateAnisotropicCoefficients(tsg::TasmanianSG; type, output)
    if getNumOutputs(tsg) == 0
        throw(TasmanianInputError("ERROR: cannot set refinement for grid with outputs = 0"))
    end
    if getNumLoaded(tsg) == 0
        throw(TasmanianInputError("ERROR: cannot call estimateAnisotropicCoefficients for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    if output == -1
        if !isSequence(tsg)
            throw(TasmanianInputError("ERROR: output = -1 can be used only for sequence grids"))
        end
    end
    if output < -1
        throw(TasmanianInputError("ERROR: output should be -1 or a non-negative integer"))
    end
    if output >= getNumOutputs(tsg)
        throw(TasmanianInputError("ERROR: output cannot exceed the index of the last output $(getNumOutputs(tsg))"))
    end
    if !(type in GlobalTypes)
        throw(TasmanianInputError("ERROR: invalid type, see TasmanianSG.GlobalTypes for list of accepted types"))
    end
    
    NumCoeffs = getNumDimensions(tsg)
    if occursin("curved", type)
        NumCoeffs = NumCoeffs * 2
    end
    
    coeff = Vector{Int32}(undef, NumCoeffs)
    tsgEstimateAnisotropicCoefficientsStatic(tsg.pGrid, type, output, coeff)
    return coeff
end
    
"""
    setSurplusRefinement!(tsg::TasmanianSG; tolerance::Float64, output::Int=-1, refinement_type::AbstractString="", level_limits=[], scale_correction = [])

using hierarchical surplusses as an error indicator, the surplus
refinement adds points to the grid to improve accuracy

when using sequence grids: this algorithm corresponds to the
                           greedy Knapsack problem

when using local polynomial or wavelet grids, this call
                         corresponds to local spatial refinement

tolerance: float (non-negative)
           the relative error tolerance, i.e.,
           we refine only for points associated with surplus
           that exceeds the tolerance

output: int (indicates the output to use)
         selects which output to use for refinement
         sequence and local polynomial grids accept -1 to
         indicate all outputs

refinement_type: hierarchical and direction refinement strategy
                 'classic'  'parents'   'direction'   'fds'   'stable'
                 applicable only for Local Polynomial and Wavelet grids

scale_correction: matrix of non-negative numbers
                  Instead of comparing the normalized surpluses to
                  the tolerance, the scaled surplus will be used.
                  The correction allows to manually guide the
                  refinement process.
                  The surplus of the j-th output of the i-th point
                  will be scaled by scale_correction[i][j].
                  scale_correction.shape[0] must be equal to
                  getNumLoaded()
                  If empty, the scale is assumed 1.0
                  scale_correction.shape[1] must be equal to the
                  number of outputs used in the process, which is
                  equal to getNumOutputs() for output == -1,
                  or 1 if output > -1.
"""
function setSurplusRefinement!(tsg::TasmanianSG; tolerance::Float64, output::Int=-1, refinement_type::AbstractString="", level_limits = [], scale_correction = [])
    if (isGlobal(tsg))
        if !(getRule(tsg) in SequenceRules)
            throw(TasmanianInputError("ERROR: setSurplusRefinement cannot be used with global grids with non-sequence rule"))
        end
    end
    if getNumLoaded(tsg) == 0
        throw(TasmanianInputError("ERROR: cannot call setSurplusRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    if tolerance < 0
        throw(TasmanianInputError("ERROR: tolerance needs to be a non-negative number"))
    end

    plevel_limits = check_level_limits_(level_limits, getNumDimensions(tsg))

    scale_correction_ = scale_correction
    if isempty(scale_correction)
        scale_correction_ = C_NULL
    else
        if eltype(scale_correction) != Float64
            scale_correction_ = convert(Vector{Float64}, scale_correction)
        end
        activeoutputs = getNumOutputs(tsg)
        if output > -1
            activeoutputs = 1
        end
        if ndims(scale_correction_) != 2
            throw(TasmanianInputError("ERROR: scale_correction must be a matrix, instead it has $(ndims(scale_correction_)) dimensions"))
        end
        if size(scale_correction_, 2) != getNumLoaded(tsg)
            throw(TasmanianInputError("ERROR: second dimension of scale_correction is $(size(scale_correction_, 2)) but the number of current points is $(getNumLoaded(tsg))"))
        end
        if size(scale_correction_, 1)  != activeoutputs
            throw(TasmanianInputError("ERROR: leading dimension of scale_correction is $(size(scale_correction_, 1)) but the refinement is set to use $(activeoutputs)"))
        end
    end

    if isempty(refinement_type)
        if !isSequence(tsg) && !isGlobal(tsg)
            throw(TasmanianInputError("ERROR: refinement_type must be specified"))
        else
            tsgSetGlobalSurplusRefinement(tsg.pGrid, tolerance, output, plevel_limits)
        end
    else
        if isSequence(tsg)
            throw(TasmanianInputError("ERROR: refinement_type not used for Sequence Grids"))
        elseif !(refinement_type in RefineTypes)
            throw(TasmanianInputError("ERROR: invalid refinement_type, see TasmanianSG.RefineTypes for the list of accepted types"))
        else
            tsgSetLocalSurplusRefinement(tsg.pGrid, tolerance, refinement_type, output, plevel_limits, scale_correction_)
        end
    end
end

"""
    getSurplusRefinement!(tsg::TasmanianSG; tolerance, output, criteria = "", level_limits = [], scale_correction = [])

Calls setSurplusRefinement!() on the inputs and then getNeededPoints().
"""
function     getSurplusRefinement!(tsg::TasmanianSG; tolerance, output, criteria = "", level_limits = [], scale_correction = [])
    setSurplusRefinement!(tsg, tolerance, output, criteria=criteria, level_limits=level_limits, scale_correction=scale_correction)
    return getNeededPoints(tsg)
end

"""
    clearRefinement!(tsg::TasmanianSG)

clear the last call to set***Refinement,
only works if called before the points are loaded, i.e.,
before loadNeededPoints()

if getNumNeeded() == 0, this call will have no effect
"""
clearRefinement!(tsg::TasmanianSG) = tsgClearRefinement(tsg.pGrid)

"""
    mergeRefinement!(tsg::TasmanianSG)

combines the loaded and needed points into a single grid
it also invalidates any currently loaded values, i.e., the
grid cannot be used for internal integration or interpolation
until loadNeededPoints() or setHierarchicalCoefficients()
is called (even if those have been called before)

if getNumNeeded() == 0, this call will have no effect
"""
mergeRefinement!(tsg::TasmanianSG) = tsgMergeRefinement(tsg.pGrid)

"""
    removePointsByHierarchicalCoefficient!(tsg::TasmanianSG; tolerance, output = -1, scale_correction = [], NumKeep = -1)

removes any points in the grid with relative surplus that
exceeds the tolerance or keeps the set number of points
with largest surplus

tolerance: float (positive)
           the relative surplus tolerance, i.e.,
           we keep only for points associated with surplus
           that exceeds the tolerance
           if NumKeep is positive, then tolerance is ignored

output: int (indicates the output to use)
        selects which output to consider
        accept -1 to indicate all outputs

scale_correction: vector or matrix
                  if output = -1 and getNumOutputs() > 1,
                  then using matrix with shape
                  getNumOutputs() X getNumLoaded() with
                  one weight per hierarchical coefficient
                  if output > -1, then using a vector with
                  one weight per point

NumKeep: int (positive or equal to -1)
         indicates the number of points to keep
         if set to -1 then tolerance is used as a cutoff
         if positive then the given number of points will be kept
"""
function removePointsByHierarchicalCoefficient!(tsg::TasmanianSG; tolerance, output = -1, scale_correction = [], num_new_points = -1)
    if !isLocalPolynomial(tsg)
        throw(TasmanianInputError("ERROR: calling removePointsByHierarchicalCoefficient for a grid that isn't local polynomial"))
    end
    if num_new_points == -1 && tolerance <= 0.0
        throw(TasmanianInputError("ERROR: tolerance must be a positive integer"))
    end
    if output < -1
        throw(TasmanianInputError("ERROR: output should be -1 or a non-negative integer"))
    end
    if output >= getNumOutputs(tsg)
        throw(TasmanianInputError("ERROR: output cannot exceed the index of the last output $(getNumOutputs(tsg) - 1))"))
    end
    if getNumLoaded(tsg) == 0
        throw(TasmanianInputError("ERROR: calling removePointsByHierarchicalCoefficient when no points are loaded"))
    end
    if num_new_points == 0 || num_new_points < -1 || num_new_points > getNumLoaded(tsg)
        throw(TasmanianInputError("ERROR: num_new_points should be either -1 or positive without exceeding the number of loaded points."))
    end
    
    if !isempty(scale_correction)
        if ndims(scale_correction) == 1
            ncol = 1
            nrow = length(scale_correction)
        elseif ndims(scale_correction) == 2
            nrow, ncol = size(scale_correction)
        else
            throw(TasmanianInputError("ERROR: scale_correction should be a vector or a matrix"))
        end
        if ncol != getNumLoaded(tsg)
            throw(TasmanianInputError("ERROR: size(scale_correction, 2) should match getNumLoaded()"))
        end
        if output == -1
            if ndims(scale_correction) != 2 && getNumOutputs(tsg) > 1
                throw(TasmanianInputError("ERROR: if output == -1 and getNumOuputs() > 1, scale_correction should be a matrix"))
            end
            if size(scale_correction, 1) != getNumOutputs(tsg)
                throw(TasmanianInputError("ERROR: size(scale_correction, 1) should be equal to getNumOutputs()"))
            end
        elseif ndims(scale_correction) != 1 && length(scale_correction) != getNumLoaded(tsg) * length(output)
            throw(TasmanianInputError("ERROR: scale_correction should be a vector of size getNumOutputs * length(output)"))
        end
    end
    
    if isempty(scale_correction)
        if num_new_points == -1
            tsgRemovePointsByHierarchicalCoefficient(tsg.pGrid, tolerance, output, C_NULL)
        else
            tsgRemovePointsByHierarchicalCoefficientHardCutoff(tsg.pGrid, num_new_points, output, C_NULL)
        end
    else
        scale_correction_ = scale_correction
        if eltype(scale_correction) != Float64
            scale_correction_ = convert(Vector{Float64}, scale_correction)
        end
        if num_new_points == -1
            tsgRemovePointsByHierarchicalCoefficient(tsg.pGrid, tolerance, output, scale_correction_)
        else
            tsgRemovePointsByHierarchicalCoefficientHardCutoff(tsg.pGrid, num_new_points, output, scale_correction_)
        end
    end
end

"""
    getHierarchicalCoefficients(tsg::TasmanianSG)

For global grids, this just returns the values loaded using the
call to loadNeededPoints().
In all other cases, this returns the list of hierarchical
coefficients, i.e., surpluses.

returns a matrix getNumOutputs() by getNumPoints()
"""
function getHierarchicalCoefficients(tsg::TasmanianSG)::Union{Matrix{Float64}, Matrix{ComplexF64}}
    NumOuts = getNumOutputs(tsg)
    if NumOuts == 0
        return zeros(0, 0)
    end
    NumPoints = getNumLoaded(tsg)
    if (NumPoints == 0)
        return zeros(NumOuts, 0)
    end
    if !isFourier(tsg)
        surplus = Matrix{Float64}(undef, NumOuts, NumPoints)
        tsgGetHierarchicalCoefficientsStatic(tsg.pGrid, surplus)
        return surplus
    else
        surplus = Matrix{Float64}(undef, NumOuts, 2 * NumPoints)
        tsgGetHierarchicalCoefficientsStatic(tsg.pGrid, surplus)
        @views return complex.(surplus[:, 1:NumPoints], surplus[:, NumPoints .+ (1:NumPoints)])
    end
end

"""
    evaluateHierarchicalFunctions(tsg::TasmanianSG, x)

evaluates the hierarchical functions at a set of points in the
domain and return a matrix with the result

x: a vector of length tsg.dimension or a matrix of dimensions tsg.dimension by getNumPoints(tsg)
    the columns indicate the points for evaluating the weights

output: returns a matrix of dimensions 
        shape == [size(x, 1), getNumPoints(tsg)]
        the values of the basis functions at the points
"""
function evaluateHierarchicalFunctions(tsg::TasmanianSG, x)
    if ndims(x) > 2
        throw(TasmanianInputError("ERROR: calling evaluateHierarchicalFunctions x should be a matrix or a vector"))
    end
    if size(x, 1) != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: calling evaluateHierarchicalFunctions size(x, 1) is not equal to getNumDimensions()"))
    end
    NumX = size(x, 2)
    if !isFourier(tsg)
        Result = Matrix{Float64}(undef, getNumPoints(tsg), NumX)
        tsgEvaluateHierarchicalFunctions(tsg.pGrid, x, NumX, Result)
        return Result
    else
        Result = Matrix{Float64}(undef, 2 * getNumPoints(tsg), NumX)
        tsgEvaluateHierarchicalFunctions(tsg.pGrid, x, NumX, Result)
        @views return complex.(Result[1:2:end, :], Result[2:2:end, :])
    end
end

"""
    getHierarchicalSupport(tsg::TasmanianSG)

returns the support of the hierarchical basis in a matrix

- only local-polynomial and wavelet grids have restricted support
- the support of all basis is restricted to the domain even
  if the support includes additional ares

output: a matrix of dimensions getNumDimensions() by getNumPoints()
        the support entries correspond to the output of getPoints()
        if x represents a point in the domain and
        if abs( x[j] - getPoints()[j, i] ) > getSupport()[j, i]
        then the i-th entry in evaluateHierarchicalFunctions(x)
        corresponding to x is guaranteed to be zero
"""
function getHierarchicalSupport(tsg::TasmanianSG)
    if getNumPoints(tsg) == 0
        return zeros(0, 0)
    end
    Result = Matrix{Float64}(undef, getNumDimensions(tsg), getNumPoints(tsg))
    tsgGetHierarchicalSupportStatic(tsg.pGrid, Result)
    return Result
end

"""
    evaluateSparseHierarchicalFunctions(tsg::TasmanianSG, x)

evaluates the hierarchical functions at a set of points in the
domain. The distinction between this function and
evaluateHierarchicalFunctions() lies in the type of the returned
result, namely a sparse vs a dense matrix.

The motivation for this function is that Local Polynomial and
Wavelet grids usually result in sparse matrices

x: a matrix with dimensions getNumDimensions() by NumX
      the entries indicate the points for evaluating

output: returns a TasmanianSimpleSparseMatrix class
        which is a simple class with three fields:
        aPntr, aIndx, and aVals which are numpy.ndarray of types
        int32, int32, and float64
        NumRows and NumCols are meta fields and have values
        NumRows = x.shape[0]
        NumCols = getNumPoints()
        The sparse matrix is compressed along the x.shape[0]
        dimension, i.e., using column compressed format
"""
function     evaluateSparseHierarchicalFunctions(tsg::TasmanianSG, x)
    if ndims(x) != 2
        throw(TasmanianInputError("ERROR: calling evaluateSparseHierarchicalFunctions(), x should be a matrix"))
    end
    if size(x, 1) != getNumDimensions(tsg)
        throw(TasmanianInputError("ERROR: calling evaluateSparseHierarchicalFunctions(), size(x, 1) is not equal to getNumDimensions()"))
    end
    NumX = size(x, 2)
    NumNZ = tsgEvaluateSparseHierarchicalFunctionsGetNZ(tsg.pGrid, x, NumX)
    Pntr = Vector{Int32}(undef, NumX + 1)
    Indx = Vector{Int32}(undef, NumNZ)
    Vals = Vector{Float64}(undef, isFourier(tsg) ? 2 * NumNZ : NumNZ)
    NumCols = NumX
    NumRows = getNumPoints(tsg)
    tsgEvaluateSparseHierarchicalFunctionsStatic(tsg.pGrid, x, NumX, Pntr, Indx, Vals)
    if isFourier(tsg)
        CVals = complex.(Vals[1:2:end], Vals[2:2:end])
        return SparseMatrixCSC(NumRows, NumCols, Pntr, Indx, CVals)
    else
        return SparseMatrixCSC(NumRows, NumCols, Pntr .+ 1, Indx .+ 1, Vals)
    end
end

"""
    setHierarchicalCoefficients!(tsg::TasmanianSG, coefficients)

 Local polynomial, Wavelet, and Sequence grids construct
 approximation using hierarchical coefficients based on the
 loaded values. This function does the opposite, the hierarchical
 coefficients are loaded directly and the values are computed
 based on the coefficients. The coefficients can be computed,
 e.g., by solving least-squares or compressed sensing problem
            min || A c - f ||
 where A is a matrix returned by evaluateHierarchicalFunctions()
 or evaluateSparseHierarchicalFunctions() for a set of points
 x; f are the values of the target function at the x
 points; and c is the vector with corresponding hierarchical
 coefficients.

 If there is a pending refinement, i.e., getNumLoaded() != 0 and
 getNumNeeded() != 0, then the refinement is discarded (since it
 was computed based on the old and now obsolete values)

 coefficients: a matrix with dimensions getNumOutputs() by getNumPoints()
               each column corresponds to the values of the
               coefficients at the corresponding point.
               The order and leading dimension must match the
               points obtained form getPoints(), the same
               order as the second dimension of
               evaluateHierarchicalFunctions().
               This matrix must be of type Matrix{ComplexF64} when
               using a Fourier grid.
"""
function setHierarchicalCoefficients!(tsg::TasmanianSG, coefficients)
    if ndims(coefficients) != 2
        throw(TasmanianInputError("ERROR: coefficients should be a matrix, instead it has $(ndims(coefficients)) dimensions"))
    end
    if size(coefficients, 2) != getNumPoints(tsg)
        throw(TasmanianInputError("ERROR: number of columns of coefficients is $(size(coefficients, 2)) but the number of current points is $(getNumNeeded(tsg))"))
    end
    if size(coefficients, 1) != getNumOutputs(tsg)
        throw(TasmanianInputError("ERROR: number of rows of coefficients is $(size(coefficients, 1)) but the number of outputs is set to $(getNumOutputs(tsg))"))
    end
    if isFourier(tsg) && !isa(coefficients, Matrix{ComplexF64})
        throw(TasmanianInputError("ERROR: using Fourier grid but coefficients is not Matrix{ComplexF64}"))
    end

    NumPoints, NumDims = size(coefficients)

    if isFourier(tsg)
        coefficientsTmp = hcat(real(coefficients), imag(coefficients))
        tsgSetHierarchicalCoefficients(tsg.pGrid, coefficientsTmp)
    else
        tsgSetHierarchicalCoefficients(tsg.pGrid, coefficients)
    end
end

"""
    integrateHierarchicalFunctions(tsg::TasmanianSG)

Computes the integrals of the hierarchical basis functions,
i.e., the same functions computed by evaluateHierarchicalFunctions().

returns a vector of length getNumPoints()
"""
function integrateHierarchicalFunctions(tsg::TasmanianSG)
    NumPoints = getNumPoints(tsg)
    if NumPoints == 0
        return zeros(0)
    end
    integrals = Vector{Float64}(undef, NumPoints)
    tsgIntegrateHierarchicalFunctionsStatic(tsg.pGrid, integrals)

    return integrals
end

"""
    getGlobalPolynomialSpace(tsg::TasmanianSG, interpolation)

returns a matrix corresponding to the polynomial space that is
integrated or interpolated exactly by the current grid

interpolation: boolean
        indicates whether to give the space associated
        with integration or interpolation

output: is a matrix of integers
        size(output, 1) is equal to iDimension
        size(output, 2) indicates the cardinality of the space
        each columns corresponds to a multi-index associated with a
        polynomial in a hierarchical tensor basis, for example,
        monomials

        see the manual for details
"""
function getGlobalPolynomialSpace(tsg::TasmanianSG, interpolation)
    interpolation_ = interploation ? 1 : 0
    NumIndexes = Ref{Cint}()
    indexes = tsgPythonGetGlobalPolynomialSpace(tsg.pGrid, interpolation_, NumIndexes)
    NumDimensions = getNumDimensions(tsg)
    polynomials = Matrix{Int}(undef, NumDimensions, NumIndexes[])
    for iI in 1:NumIndexes[]
        for iJ in 1:NumDimensions
            polynomials[iJ, iI] = indexes[(iI - 1)*NumDimensions + iJ]
        end
    end
    tsgDeleteInts(indexes)
    return polynomials
end

"""
    enableAcceleration!(tsg::TasmanianSG, acceleration_type, GPUID = 0)

Enables the use of accelerated backend libraries and extensions,
such as BLAS and CUDA.
Each acceleration type requires corresponding CMake compile
options, otherwise the backend will fallback to the closest
available options.

sAccelerationType: string

  'none'
      core fallback mode, relies on sequential implementation
      if compiled with Tasmanian_ENABLE_OPENMP this will use
      simple "omp parallel for" to take advantage of multiple cpu cores

  'cpu-blas'
      uses BLAS level 2 and 3 functions for acceleration of batch
      evaluations
      requires Tasmanian_ENABLE_BLAS=ON
      this is the default mode, if available

  'gpu-default'
      uses CUDA kernels, cuBlas, cuSparse and MAGMA libraries for
      accelerated matrix operations, e.g., cublasDgemm
      refer to TasGrid::TypeAcceleration for more details

  'gpu_cublas'
      uses the Nvidia cuBlas and cuSparse libraries

  'gpu-cuda'
      uses custom CUDA kernels in addition to the accelerated
      linear algebra libraries

   'gpu-magma'
      uses the custom CUDA kernels and the MAGMA library in place
      of the default Nvidia libraries

GPUID: integer
      indicates the GPU device to use, if set to None then device
      zero will be used first or the device set with setGPUID()
"""
function enableAcceleration!(tsg::TasmanianSG, acceleration_type; GPUID::Union{Int, Nothing} = nothing)
    if !(acceleration_type in AccelTypes)
        throw(TasmanianInputError("ERROR: invalid acceleration type"))
    end
    if isnothing(GPUID)
        tsgEnableAcceleration(tsg.pGrid, acceleration_type)
    else
        if GPUID < 0 || GPUID >= getNumGPUs(tsg)
            throw(TasmanianInputError("ERROR: invalid GPU ID number"))
        end
        tsgEnableAcceleration(tsg.pGrid, acceleration_type, GPUID)
    end
end

"""
    getAccelerationType(tsg::TasmanianSG)
returns the type of acceleration set by enableAcceleration
"""
function getAccelerationType(tsg::TasmanianSG)
    acceleration_type = tsgGetAccelerationType(tsg.pGrid)
    iend = 0
    while unsafe_load(acceleration_type + iend) > 0
        iend += 1
    end
    return unsafe_string(acceleration_type)
end

"""
    isAccelerationAvailable(tsg::TasmanianSG, acceleration_type)

returns true if the library has been compiled with support
for acceleration_type.
Even if this returns false, you can use the type for
enableAcceleration, but the library will default to the next
available type (see the Manual)
"""
function isAccelerationAvailable(tsg::TasmanianSG, acceleration_type)
    if !(acceleration_type in AccelTypes)
            throw(TasmanianInputError("ERROR: invalid acceleration type"))
    end
    return convert(Bool, tsgIsAccelerationAvailable(acceleration_type))
end

"""
    setGPUID!(tsg::TasmanianSG, GPUID)

when using cuda on a machine with multiple GPUs, this helps set
the GPU for this grid
NOTE: each instance of the sparse grids class holds a separate
      instance of GPUID and different grids can be assigned to
      different GPUs (on multi-gpu system)
GPUID can be changed at any time, however, this will cause
some of the internal cache to be invalidated and it may lead
to extraneous data movement

calling read or make***Grid will reset the selected GPU

defaults to 0

this doesn't do anything unless enableAcceleration is called
using a "gpu-" acceleration type
"""
function setGPUID!(tsg::TasmanianSG, GPUID)
    if GPUID < 0 || GPUID >= getNumGPUs()
        throw(TasmanianInputError("ERROR: invalid GPU ID number"))
    end
    tsgSetGPUID(tsg.pGrid, GPUID)
end

"""
    getGPUID(tsg::TasmanianSG)

returns the GPU ID set using setGPUID
"""
getGPUID(tsg::TasmanianSG) = tsgGetGPUID(tsg.pGrid)

"""
    getNumGPUs()
        returns the number of available GPUs according to cuda

        this is one of several functions designed to allow basic
        management of multi-gpu setup with only Tasmanian module
"""
getNumGPUs() = tsgGetNumGPUs()

"""
    getGPUMemory(tsg::TasmanianSG, GPUID)

returns the total memory (in MegaBytes, 1024**2 bytes) of the
corresponding GPU

this is one of several functions designed to allow basic
management of multi-gpu setup with only Tasmanian module
"""
function getGPUMemory(tsg::TasmanianSG, GPUID)
    if GPUID < 0 || GPUID >= getNumGPUs()
        throw(TasmanianInputError("ERROR: invalid GPU ID number"))
    end
    return tsgGetGPUMemory(GPUID)
end

"""
    getGPUName(tsg::TasmanianSG, GPUID)

return the cuda name ID of the corresponding GPU

this is one of several functions designed to allow basic
management of multi-gpu setup with only Tasmanian module
"""
function getGPUName(tsg::TasmanianSG, GPUID)
    if GPUID < 0 || GPUID >= getNumGPUs()
        throw(TasmanianInputError("ERROR: invalid GPU ID number"))
    end
    buffer_size = 256
    name = Vector{Cchar}(undef, buffer_size)
    NumChars = Ref{Int32}()
    tsgGetGPUName(GPUID, buffer_size, name, NumChars)
    return String([Char(n) for n in name if n != 0])
end

"""
    printStats(tsg::TasmanianSG)

calls the library printStats() function, which displays basic
information about this instance of the grid
"""
printStats(tsg::TasmanianSG) = tsgPrintStats(tsg.pGrid)

#=
"""
    plotPoints2D(tsg::TasmanianSG, pAxisObject=tsgPlot, sStyle="bo", iMarkerSize=3)

plots the points in a 2D plot using matplotlib.pyplot
applicable only for grids with iDimensions == 2

pAxisObject: axis object from the matplotlib.pyplot package

sStyle: string
        the matplotlib.pyplot style, e.g.,
        'ko' will make black cirlces, 'rx' will use red crosses

iMarkerSize: positive integer
             the marker size for plotting the points
"""
function plotPoints2D(tsg::TasmanianSG, pAxisObject=tsgPlot, sStyle="bo", iMarkerSize=3)
    if getNumDimensions(tsg) != 2
        throw(TasmanianInputError("ERROR: cannot plot a grid with other than 2 dimensions"))
    end

    aPoints = getPoints(tsg)

    plot(aPoints[:, 1], aPoints[:, 2])
end

"""
    plotResponse2D(tsg::TasmanianSG, output=0, iNumDim0=100, iNumDim1=100, pAxisObject=tsgPlot, sCmap="jet")

plots the response in a 2D plot using matplotlib.pyplot
applicable only for grids with iDimensions == 2

output is the output to use for plotting

iNumDim0, iNumDim1: positive integers
       the points for the plot are selected on a dense grid with
       number of points iNumDim0 and iNumDim1 in dimensions
       0 and 1 respectively

pAxisObject: axis object from the matplotlib.pyplot package

sCmap: string indicating the map to use, e.g., "jet" or "heat"
"""
function plotResponse2D(tsg::TasmanianSG, output=0, iNumDim0=100, iNumDim1=100, pAxisObject=tsgPlot, sCmap="jet")
        if getNumDimensions(tsg) != 2
            throw(TasmanianInputError("ERROR: cannot plot a grid with other than 2 dimensions"))
        end
        if (iNumDim0 < 1)
            throw(TasmanianInputError("ERROR: the number of points in dimension 1 should be at least 1"))
        end
        if (iNumDim1 < 1)
            throw(TasmanianInputError("ERROR: the number of points in dimension 2 should be at least 1"))
        end     
        aPoints = getPoints()

        xmin = min(aPoints[:,0])
        xmax = max(aPoints[:,0])
        fymin = min(aPoints[:,1])
        fymax = max(aPoints[:,1])
        if (xmin == xmax)
            xmin = xmin - 0.1
            xmax = xmax + 0.1
        if (fymin == fymax)
            fymin = fymin - 0.1
            fymax = fymax + 0.1

        x = np.linspace(xmin, xmax, iNumDim0)
        y = np.linspace(fymax, fymin, iNumDim1) # flip the order of y to match the top-to-bottom pixel indexing

        XX, YY = np.meshgrid(x, y)
        ZZ = evaluateBatch(np.vstack((XX.reshape((iNumDim0*iNumDim1,)), YY.reshape((iNumDim0*iNumDim1,)))).T)
        ZZ = ZZ[:,output].reshape((iNumDim0,iNumDim1))

        pAxisObject.imshow(ZZ, cmap=sCmap, extent=[xmin, xmax, fymin, fymax])
        end
=#                  




    
