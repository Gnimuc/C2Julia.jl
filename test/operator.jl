using C2Julia
using Clang
using Test

@testset "unary operator" begin
    tu = parse_header(joinpath(@__DIR__, "c", "operator.h"))
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        child_cursors = children(root_cursor)
        op_func = search(root_cursor, "operators")[1]
        body = children(op_func)[1]
        op = children(body)[4]
        @test string(translate(op)) == "+5"
        op = children(body)[5]
        @test string(translate(op)) == "-5"
        op = children(body)[6]
        @test translate(op).expr == :(!z)
        op = children(body)[7]
        @test translate(op).expr == :(Cint(z))
        op = children(body)[8]
        @test string(translate(op)) == "@+ y"
        op = children(body)[9]
        @test translate(op).expr == :(y+=1)
        op = children(body)[10]
        @test string(translate(op)) == "@- x"
        op = children(body)[11]
        @test translate(op).expr == :(x-=1)
        op = children(body)[12]
        @test translate(op)[].expr == :(xp = Ref(x))
        op = children(body)[13]
        @test string(translate(op)) == "xp[] = 1"
        op = children(body)[14]
        @test translate(op).expr == :(Bool(z+1))
    end
end
