using C2Julia
using Clang
using Test

@testset "array" begin
    tu = parse_header(joinpath(@__DIR__, "c", "array.h"))
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        func = search(root_cursor, "array")[1]
        body = children(func)[1]
        arrayinit = children(body)[1]
        @test string(translate(arrayinit)[]) == "a = [0, 1, 2, 3, 4, 5, 6, 7, 8]"
        arrayinit2 = children(body)[2]
        @test string(translate(arrayinit2)[]) == "b = [[0, 1, 2, 3], [4, 5, 6, 7]]"
        indexing = children(body)[3]
        @test string(translate(indexing)) == "(b[1])[2] = a[3]"
    end
end
