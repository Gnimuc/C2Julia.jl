translate(cursor::CLDeclRefExpr) = translate(getref(cursor))
translate(cursor::CLVarDecl) = Meta.parse(spelling(cursor))

function translate(cursor::CLUnexposedExpr)
    child_cursors = children(cursor)
    if length(child_cursors) == 1
        translate(child_cursors[1])
    else
        @warn "translate subroutine for $cursor is not implemented." dumpobj(cursor)
        return Expr(:null)
    end
end
