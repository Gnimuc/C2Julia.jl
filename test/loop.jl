using C2Julia
using Clang
using Test

@testset "loops" begin
    tu = parse_header(joinpath(@__DIR__, "c", "loop.h"), includes=[LLVM_INCLUDE])
    @testset "while-loop" begin
        GC.@preserve tu begin
            root_cursor = getcursor(tu)
            loop_func = search(root_cursor, "loop_stmt")[1]
            body = children(loop_func)[1]
            while_stmt = children(body)[6]
            expr = :(while k < 3
                        x = 0
                        k = k + 1
                    end)
            x = Base.remove_linenums!(expr)
            y = translate(while_stmt)
            @test string(x) == string(y)
        end
    end

    @testset "do-while-loop" begin
        GC.@preserve tu begin
            root_cursor = getcursor(tu)
            loop_func = search(root_cursor, "loop_stmt")[1]
            body = children(loop_func)[1]
            do_stmt = children(body)[8]
            expr = :(while true
                        m = m + 1
                        m < 7 || break
                    end)
            x = Base.remove_linenums!(expr)
            y = translate(do_stmt)
            @test string(x) == string(y)
        end
    end

    @testset "for-loop" begin
        GC.@preserve tu begin
            root_cursor = getcursor(tu)
            loop_func = search(root_cursor, "loop_stmt")[1]
            body = children(loop_func)[1]
            for_stmt = children(body)[4]
            expr = Expr(:macrocall, Symbol("@cfor"), nothing,
                        :(j=0), :(j<10), :(j=j+1), Expr(:block, :(x = x - 1)))
            x = Base.remove_linenums!(expr)
            y = translate(for_stmt)
            @test string(x) == string(y)

            infinite_stmt = children(body)[10]
            expr = Expr(:macrocall, Symbol("@cfor"), nothing,
                        nothing, nothing, nothing, Expr(:block, :(x = 0)))
            x = Base.remove_linenums!(expr)
            y = translate(infinite_stmt)
            @test string(x) == string(y)
        end
    end
end
