translate(cursor::CLDeclRefExpr) = MetaExpr(Symbol(spelling(cursor)), cursor)

function getranges(cursor::CLCursor)
    tu = clang_Cursor_getTranslationUnit(cursor)
    source_range = clang_getCursorExtent(cursor)
    range_start = clang_getRangeStart(source_range)
    range_end = clang_getRangeEnd(source_range)

    file = Ref{CXFile}(C_NULL)
    start_line = Ref{Cuint}(0)
    end_line = Ref{Cuint}(0)
    start_column = Ref{Cuint}(0)
    end_column = Ref{Cuint}(0)
    offset = Ref{Cuint}(0)
    clang_getExpansionLocation(range_start, file, start_line, start_column, offset)
    clang_getExpansionLocation(range_end, file, end_line, end_column, offset)

    expanded_start = clang_getLocation(tu, file[], start_line[], start_column[])
    expanded_end = clang_getLocation(tu, file[], end_line[], end_column[])
    expanded_range = clang_getRange(expanded_start, expanded_end)

    return source_range, expanded_range
end

function translate(cursor::CLParenExpr)
    source, expanded = getranges(cursor)
    if source == expanded
        return translate(first(children(cursor)))
    else
        # we do not translate function-like macros
        tu = clang_Cursor_getTranslationUnit(cursor)
        toks = Clang.TokenList(tu, expanded)
        code = mapreduce(x->x.text, *, collect(toks))
        expr = Expr(:macrocall, Symbol("@error"), nothing, code)
        return MetaExpr(expr, cursor)
    end
end

function translate(cursor::CLMemberRefExpr)
    sym = Symbol(spelling(cursor))
    meta = translate(first(children(cursor)))
    meta.expr = Expr(:(.), meta.expr, Meta.quot(sym))
    return meta
end

function translate(cursor::CLCStyleCastExpr)
    child_cursors = children(cursor)
    jltype_sym = clang2julia(type(cursor))
    if jltype_sym == :Cvoid
        return translate(last(child_cursors))  # don't cast anything to Nothing
    else
        meta = translate(last(child_cursors))
        meta.expr = Expr(:call, jltype_sym, meta.expr)
        return meta
    end
end

function translate(cursor::CLUnexposedExpr)
    child_cursors = children(cursor)
    # if length(child_cursors) == 1
    return translate(first(child_cursors))
    # else
    #     @warn "translate subroutine for $cursor is not implemented."
    #     return Expr(:null)
    # end
end

function translate(cursor::CLUnaryExpr)
    child_cursors = children(cursor)
    # if length(child_cursors) == 1
    return translate(first(child_cursors))
    # else
    #     @warn "translate subroutine for $cursor is not implemented."
    #     return Expr(:null)
    # end
end

function translate(cursor::CLCallExpr)
    child_cursors = children(cursor)
    func_name = spelling(cursor)
    func_expr = Expr(:call, Symbol(func_name))
    meta = MetaExpr[]
    for c in child_cursors[2:end]
        meta = translate(c)
        push!(func_expr.args, meta.expr)
    end
    return MetaExpr(func_expr)
end
