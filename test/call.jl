using C2Julia
using C2Julia.CSyntax
using Clang
using Test

@testset "call" begin
    tu = parse_header(joinpath(@__DIR__, "c", "call.h"), includes=[LLVM_INCLUDE])
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        caller_cursor = search(root_cursor, "caller")[]
        callee_cursor = search(root_cursor, "callee")[]

        translate(caller_cursor) |> println
        translate(callee_cursor) |> println

        string(translate(caller_cursor)) |> Meta.parse |> eval
        string(translate(callee_cursor)) |> Meta.parse |> eval

        @test caller() == 0
    end
end
