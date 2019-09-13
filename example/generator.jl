using C2Julia
using Clang
using Clang.LibClang

liblsoda_include = joinpath(@__DIR__, "liblsoda", "src")
liblsoda_c = joinpath(@__DIR__, "liblsoda", "src", "lsoda.c")

ctx = DefaultContext()

std_headers = find_std_headers()
compiler_flags = ["-I$liblsoda_include", "-I$LLVM_INCLUDE", ("-I".*std_headers)...]

push!(ctx.trans_units, TranslationUnit(ctx.index,
                                       liblsoda_c,
                                       compiler_flags,
                                       CXTranslationUnit_None))

root_cursor = getcursor(ctx.trans_units[1])
header = spelling(root_cursor)
ctx.children = children(root_cursor)

check_opt = search(root_cursor, "check_opt")[1]

outpath = joinpath(@__DIR__, "lsoda_test.jl")
open(outpath, "w+") do f
    println(f, translate(check_opt).expr)
end
