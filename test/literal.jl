using C2Julia
using Clang
using Test

@testset "literal trivial" begin
    tu = parse_header(joinpath(@__DIR__, "c", "literal.h"))
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        child_cursors = children(root_cursor)
        decl = child_cursors[5]
        integer = children(decl)[]
        @test translate(integer).expr == 0

        decl = child_cursors[6]
        float = children(decl)[]
        @test translate(float).expr == :(Float32(1.0))

        decl = child_cursors[7]
        hex = children(decl)[]
        @test translate(hex).expr == :(0x10)

        decl = child_cursors[8]
        hexu = children(decl)[]
        @test translate(hexu).expr == :(UInt32(0x10))

        decl = child_cursors[9]
        oct = children(decl)[]
        @test translate(oct).expr == :(0o123)

        decl = child_cursors[10]
        char = children(decl)[]
        @test translate(char).expr == :('f')

        decl = child_cursors[11]
        str = children(decl)[]
        @test translate(str).expr == :("abcdefg")
    end
end
