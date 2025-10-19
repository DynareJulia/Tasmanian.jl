using Clang
using Clang.Generators
using Clang.LibClang.Clang_jll
using Tasmanian_jll

cd(@__DIR__)

include_dir = normpath(Tasmanian_jll.artifact_dir, "include")

# wrapper generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# add compiler flags, e.g. "-DXXXXXXXXX"
args = get_default_args()
push!(args, "-I$include_dir")

# only wrap libclang headers in include/clang-c
header_dir = include_dir
@show include_dir
headers = [joinpath(header_dir, header)
           for header in readdir(header_dir) if endswith(header, ".h")]
@show headers

# create context
ctx = create_context(headers, args, options)

# run generator
build!(ctx)
