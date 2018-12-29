module C2Julia

using Clang
using Clang.LibClang
import Clang: TokenList

translate(cursor::CLCursor) = "not implemented yet"

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
