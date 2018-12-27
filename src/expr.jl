translate(cursor::CLDeclRefExpr) = Symbol(spelling(cursor))
translate(cursor::CLParenExpr) = translate(children(cursor)[])

function translate(cursor::CLMemberRefExpr)
    sym = Symbol(spelling(cursor))
    return Expr(:(.), translate(children(cursor)[]), Meta.quot(sym))
end

function translate(cursor::CLCStyleCastExpr)
    child_cursors = children(cursor)
    jltype_sym = clang2julia(type(cursor))
    if jltype_sym == :Cvoid
        return translate(last(child_cursors))  # don't cast anything to Nothing
    else
        return Expr(:call, jltype_sym, translate(last(child_cursors)))
    end
end

function translate(cursor::CLUnexposedExpr)
    child_cursors = children(cursor)
    if length(child_cursors) == 1
        translate(child_cursors[])
    else
        @warn "translate subroutine for $cursor is not implemented." dumpobj(cursor)
        return Expr(:null)
    end
end

function translate(cursor::CLUnaryExpr)
    child_cursors = children(cursor)
    if length(child_cursors) == 1
        translate(child_cursors[])
    else
        @warn "translate subroutine for $cursor is not implemented." dumpobj(cursor)
        return Expr(:null)
    end
end

function translate(cursor::CLCallExpr)
    child_cursors = children(cursor)
    func_name = spelling(cursor)
    func_expr = Expr(:call, Symbol(func_name))
    for c in child_cursors[2:end]
        push!(func_expr.args, translate(c))
    end
    return func_expr
end
