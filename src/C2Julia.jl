module C2Julia

using Clang
using Clang.LibClang
import Clang: TokenList

translate(cursor::CLCursor) = "not implemented yet"

include("literal.jl")
include("operator.jl")
include("expr.jl")
include("stmt.jl")

export translate

end # module
