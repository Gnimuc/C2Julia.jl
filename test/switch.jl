using C2Julia
using CSyntax
using Clang
using Test

@testset "switch" begin
    tu = parse_header(joinpath(@__DIR__, "c", "switch.h"), includes=[LLVM_INCLUDE])
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        switch_func = search(root_cursor, "switch_func")[1]
        body = children(switch_func)[1]
        switch = children(body)[3]
        x = translate(switch)
        expr = Expr(:macrocall, Symbol("@switch"), nothing, :str, Expr(:block,
                    Expr(:macrocall, Symbol("@case"), nothing, 'a'),
                    :(x = 1),
                    :(x = 0),
                    Expr(:break),
                    Expr(:macrocall, Symbol("@case"), nothing, 'b'),
                    Expr(:if, :(x == 2), :(x = 7)),
                    Expr(:macrocall, Symbol("@case"), nothing, 'c'),
                    Expr(:macrocall, Symbol("@default"), nothing),
                    :(x = 10)
                    ))
        y = Base.remove_linenums!(expr)
        @test string(x) == string(y)
   end
end
