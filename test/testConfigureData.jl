artifact_dir = Tasmanian_jll.artifact_dir
TasmanianConfig_hpp = read(joinpath(artifact_dir, "include",
                                    "TasmanianConfig.hpp"), String)

bEnableSyncTests = true

sLibPath = dirname(Tasmanian_jll.get_libtasmaniansparsegrid_path())

iGPUID = -1

bHasBlas = occursin("#define Tasmanian_ENABLE_BLAS", TasmanianConfig_hpp)
bHasSycl = occursin("#define Tasmanian_ENABLE_DPCPP", TasmanianConfig_hpp)
bHasCuBlas = (bHasSycl ||
              occursin("#define Tasmanian_ENABLE_CUDA", TasmanianConfig_hpp) ||
              occursin("#define Tasmanian_ENABLE_HIP", TasmanianConfig_hpp))

bHasCuda = bHasCuBlas


bUsingMSVC = false #("@CMAKE_CXX_COMPILER_ID@" == "MSVC")

sGaussPattersonTableFile = joinpath(artifact_dir, "share",
                                    "Tasmanian",
                                    "GaussPattersonTableFile")

