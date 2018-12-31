function translate(cursor::CLVarDecl)
    child_cursors = children(cursor)
    cursor_sym = Symbol(spelling(cursor))
    if isempty(child_cursors)
        return MetaExpr(cursor_sym, cursor)
    elseif length(child_cursors) == 1
        var_meta = translate(first(child_cursors))
        var_meta.expr = Expr(:(=), cursor_sym, var_meta.expr)
        return var_meta
    else
        var_meta = translate(first(child_cursors))
        list_meta = translate(last(child_cursors))
        if Meta.isexpr(var_meta.expr, :call)
            push!(var_meta.expr.args, list_meta.expr.args...)
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
        @warn "not implemented yet"
        return MetaExpr(Expr(:null))
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
