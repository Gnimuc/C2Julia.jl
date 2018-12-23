using C2Julia
using Clang
using Test

@testset "if-then" begin
    tu = parse_header(joinpath(@__DIR__, "c", "if.h"))
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        if_func = search(root_cursor, "if_stmt")[1]
        body = children(if_func)[1]
        single_stmt = children(body)[3]
        expr = :(if 0 == a
                    b = 2
                end)
        x = Base.remove_linenums!(expr)
        y = translate(single_stmt)
        @test string(x) == string(y)
    end
end

@testset "if-else" begin
    tu = parse_header(joinpath(@__DIR__, "c", "if.h"))
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        if_func = search(root_cursor, "if_stmt")[1]
        body = children(if_func)[1]
        else_stmt = children(body)[4]
        expr = :(if 0 == a
                    b = 2
                else
                    b = 3
                end)
        x = Base.remove_linenums!(expr)
        y = translate(else_stmt)
        @test string(x) == string(y)
   end
end

@testset "if-nested" begin
    tu = parse_header(joinpath(@__DIR__, "c", "if.h"))
    GC.@preserve tu begin
        root_cursor = getcursor(tu)
        if_func = search(root_cursor, "if_stmt")[1]
        body = children(if_func)[1]
        if_elseif_stmt = children(body)[5]
        expr = :(if 0 == a
                    b = 2
                elseif a < 0
                    b = 3
                else
                    b = 4
                end)
        x = Base.remove_linenums!(expr)
        y = translate(if_elseif_stmt)
        @test string(x) == string(y)

        nested_stmt = children(body)[6]
        expr = :(if 0 == a
                    b = 2
                else
                    if a < 0
                        b = 3
                    else
                        b = 4
                    end
                end)
        x = Base.remove_linenums!(expr)
        y = translate(nested_stmt)
        @test string(x) == string(y)
    end
end
