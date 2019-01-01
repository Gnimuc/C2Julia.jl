# TODO: clean up ugly routine
function translate(cursor::CLVarDecl)
    child_cursors = children(cursor)
    cursor_sym = Symbol(spelling(cursor))
    if isempty(child_cursors)
        return MetaExpr(cursor_sym, cursor)
    elseif length(child_cursors) == 1
        var_cursor = first(child_cursors)
        var_meta = translate(var_cursor)
        var_meta.expr == Expr(:null) && return MetaExpr(cursor_sym, cursor)
        var_meta.expr = Expr(:(=), cursor_sym, var_meta.expr)
        return var_meta
    else
        var_meta = translate(first(child_cursors))
        list_meta = translate(last(child_cursors))
        if Meta.isexpr(var_meta.expr, :call) && kind(last(child_cursors)) != CXCursor_CallExpr
            push!(var_meta.expr.args, list_meta.expr.args...)
            var_meta.expr = Expr(:(=), cursor_sym, var_meta.expr)
        else
            var_meta.expr = Expr(:(=), cursor_sym, list_meta.expr)
        end
        return var_meta
    end
end

function translate(cursor::CLTypeRef)
    origin = getref(cursor)
    if kind(origin) == CXCursor_StructDecl
        return MetaExpr(Expr(:call, Symbol(spelling(origin))))
    else
        # TODO: check this out
        @warn "unexpected subroutine reached."
        return MetaExpr(Expr(:null), cursor)
    end
end

translate(cursor::CLParmDecl) = MetaExpr(Symbol(spelling(cursor)), cursor)

function translate(cursor::CLFunctionDecl)
    child_cursors = children(cursor)
    signature = Expr(:call, Symbol(spelling(cursor)))
    body = Expr(:null)
    for c in child_cursors
        meta = translate(c)
        if kind(c) == CXCursor_ParmDecl
            push!(signature.args, meta.expr)
        else
            body = meta.expr
        end
    end
    return MetaExpr(Expr(:function, signature, body))
end
