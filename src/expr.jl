translate(cursor::CLDeclRefExpr) = Symbol(spelling(cursor))
translate(cursor::CLParenExpr) = translate(children(cursor)[])

function translate(cursor::CLMemberRefExpr)
    sym = Symbol(spelling(cursor))
    return Expr(:(.), translate(children(cursor)[]), Meta.quot(sym))
end

function translate(cursor::CLCStyleCastExpr)
    child = children(cursor)[]
    jltype_sym = clang2julia(type(cursor))
    jltype_sym == :Cvoid && return translate(child)  # don't cast anything to Nothing
    Expr(:call, jltype_sym, translate(child))
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

function translate(cursor::CLCallExpr)
    child_cursors = children(cursor)
    func_name = spelling(cursor)
    func_expr = Expr(:call, Symbol(func_name))
    for c in child_cursors[2:end]
        push!(func_expr.args, translate(c))
    end
    return func_expr
end
