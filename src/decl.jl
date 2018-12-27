function translate(cursor::CLVarDecl)
    child_cursors = children(cursor)
    if isempty(child_cursors)
        return Symbol(spelling(cursor))
    else
        return Expr(:(=), Symbol(spelling(cursor)), translate(child_cursors[]))
    end
end

function translate(cursor::CLTypeRef)
    origin = getref(cursor)
    if kind(origin) == CXCursor_StructDecl
        return Expr(:call, Symbol(spelling(origin)))
    else
        @error "not implemented yet"
    end
end

translate(cursor::CLParmDecl) = Symbol(spelling(cursor))

function translate(cursor::CLFunctionDecl)
    child_cursors = children(cursor)
    signature = Expr(:call, Symbol(spelling(cursor)))
    body = Expr(:null)
    for c in child_cursors
        if kind(c) == CXCursor_ParmDecl
            push!(signature.args, translate(c))
        else
            body = translate(c)
        end
    end
    return Expr(:function, signature, body)
end
