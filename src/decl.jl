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
