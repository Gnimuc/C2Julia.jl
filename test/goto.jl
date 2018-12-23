using C2Julia
using Clang
using Test

@testset "goto" begin
    tu = parse_header(joinpath(@__DIR__, "c", "goto.h"))
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        op_func = search(root_cursor, "goto_func")[1]
        body = children(op_func)[1]
        op = children(body)[2]
        @test string(translate(op)) == "@goto jump"
        op = children(body)[3]
        @test string(translate(op)) == "@label jump"
    end
end
