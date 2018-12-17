module C2Julia

using Clang
using Clang.LibClang

include("translate.jl")
export translate

end # module
