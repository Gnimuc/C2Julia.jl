const C2JULIA_NULL_CURSOR = CLCursor(getcursor())

"""
    MetaExpr
Symbols, Exprs, Numbers, Chars, Strings
"""
mutable struct MetaExpr
    expr::Union{Symbol,Expr,Number,AbstractChar,AbstractString}
    leafcursor::CLCursor
end
MetaExpr(x) = MetaExpr(x, C2JULIA_NULL_CURSOR)

Base.length(x::MetaExpr) = 1
Base.iterate(x::MetaExpr) = (x, nothing)
Base.string(x::MetaExpr) = string(x.expr)
Base.show(io::IO, x::MetaExpr) = show(io, x.expr)
