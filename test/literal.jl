using C2Julia
using Clang
using Test

@testset "literal trivial" begin
    tu = parse_header(joinpath(@__DIR__, "c", "literal.h"), includes=[LLVM_INCLUDE])
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        child_cursors = children(root_cursor)
        decl = search(child_cursors, x->name(x)=="a")[]
        integer = children(decl)[]
        @test translate(integer).expr == 0

        decl = search(child_cursors, x->name(x)=="b")[]
        float = children(decl)[]
        @test translate(float).expr == :(Float32(1.0))

        decl = search(child_cursors, x->name(x)=="c")[]
        hex = children(decl)[]
        @test translate(hex).expr == :(0x10)

        decl = search(child_cursors, x->name(x)=="d")[]
        hexu = children(decl)[]
        @test translate(hexu).expr == :(UInt32(0x10))

        decl = search(child_cursors, x->name(x)=="e")[]
        oct = children(decl)[]
        @test translate(oct).expr == :(0o123)

        decl = search(child_cursors, x->name(x)=="f")[]
        char = children(decl)[]
        @test translate(char).expr == :('f')

        decl = search(child_cursors, x->name(x)=="g")[]
        str = children(decl)[]
        @test translate(str).expr == :("abcdefg")
    end
end
