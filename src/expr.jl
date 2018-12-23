translate(cursor::CLDeclRefExpr) = Symbol(spelling(cursor))

function translate(cursor::CLVarDecl)
    child_cursors = children(cursor)
    if isempty(child_cursors)
        return Symbol(spelling(cursor))
    else
        return Expr(:(=), Symbol(spelling(cursor)), translate(child_cursors[]))
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

function translate(cursor::CLCallExpr)
    child_cursors = children(cursor)
    func_name = spelling(cursor)
    func_expr = Expr(:call, Symbol(func_name))
    for c in child_cursors[2:end]
        push!(func_expr.args, translate(c))
    end
    return func_expr
end
