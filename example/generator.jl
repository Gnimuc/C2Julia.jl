using C2Julia
using Clang
using Clang.LibClang

lsoda = joinpath(@__DIR__, "liblsoda")
liblsoda_c = joinpath(lsoda, "lsoda.c")

ctx = DefaultContext()

std_headers = find_std_headers()
compiler_flags = ["-I$LLVM_INCLUDE", ("-I".*std_headers)...]

push!(ctx.trans_units, TranslationUnit(ctx.index,
                                       liblsoda_c,
                                       compiler_flags,
                                       CXTranslationUnit_None))

root_cursor = getcursor(ctx.trans_units[1])
header = spelling(root_cursor)
ctx.children = children(root_cursor)

funs = filter(x -> x isa CLFunctionDecl && filename(x)==header, ctx.children)
cfuns = map(x->search(root_cursor, spelling(x)) |> last, funs)

outpath = joinpath(@__DIR__, "lsoda_test.jl")
open(outpath, "w+") do f
    for cfun in cfuns
        println(f, translate(cfun).expr)
    end
end
