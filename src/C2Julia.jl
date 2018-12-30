module C2Julia

using Clang
using Clang.LibClang
import Clang: TokenList

function translate(cursor::CLCursor)
    expr = Expr(:macrocall, Symbol("@error"), nothing, "not implemented yet or unknown bugs")
    return MetaExpr(expr, cursor)
end

include("CSyntax.jl")
using .CSyntax

include("types.jl")
export MetaExpr

include("literal.jl")
include("operator.jl")
include("expr.jl")
include("stmt.jl")
include("decl.jl")

export translate

end # module
