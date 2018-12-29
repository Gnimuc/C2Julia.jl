function translate(cursor::CLVarDecl)
    child_cursors = children(cursor)
    cursor_sym = Symbol(spelling(cursor))
    if isempty(child_cursors)
        return MetaExpr(cursor_sym, cursor)
    else
        meta = translate(first(child_cursors))
        meta.expr = Expr(:(=), cursor_sym, meta.expr)
        return meta
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
