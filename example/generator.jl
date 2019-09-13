using C2Julia
using Clang
using Clang.LibClang

cvode = joinpath(homedir(), "Codes", "CVODE")
libcvode_include = joinpath(cvode, "include")
libcvode_c = joinpath(cvode, "source", "cvode.c")

ctx = DefaultContext()

std_headers = find_std_headers()
compiler_flags = ["-I$libcvode_include", "-I$LLVM_INCLUDE", ("-I".*std_headers)...]

push!(ctx.trans_units, TranslationUnit(ctx.index,
                                       libcvode_c,
                                       compiler_flags,
                                       CXTranslationUnit_None))

root_cursor = getcursor(ctx.trans_units[1])
header = spelling(root_cursor)
ctx.children = children(root_cursor)

funs = filter(x -> x isa CLFunctionDecl && filename(x)==header, ctx.children)
cfuns = map(x->search(root_cursor, spelling(x)) |> last, funs)

outpath = joinpath(@__DIR__, "cvode_test.jl")
open(outpath, "w+") do f
    for cfun in cfuns
        println(f, translate(cfun).expr)
    end
end
