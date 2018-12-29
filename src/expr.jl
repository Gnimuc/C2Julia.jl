translate(cursor::CLDeclRefExpr) = MetaExpr(Symbol(spelling(cursor)), cursor)
translate(cursor::CLParenExpr) = translate(first(children(cursor)))

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
