using C2Julia
using C2Julia.CSyntax
using Clang
using Test

include("vsSwitch.jl")
@testset "test against Switch.jl" begin
    @testset "fall through" begin
        @test switch_with_fallthrough(1) == Int[1, 2, 3, 4]
        @test switch_with_fallthrough(3) == Int[3, 4]
        @test switch_with_fallthrough(100) == Int[4]
    end
    @testset "leading default" begin
        @test switch_with_leading_default(1) == Int[1, 2, 3]
        @test switch_with_leading_default(3) == Int[3]
        @test switch_with_leading_default(100) == Int[4, 1, 2, 3]
    end
    @testset "with expressions" begin
        @test switch_with_expressions(1) == Int[1, 2, 3, 4]
        @test switch_with_expressions(3) == Int[3, 4]
        @test switch_with_expressions(100) == Int[4]
    end
    @testset "with break" begin
        @test switch_with_break(1) == Int[1]
        @test switch_with_break(3) == Int[3]
        @test switch_with_break(100) == Int[4]
    end
end

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
