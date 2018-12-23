abstract type UnaryOperator end

for op in [:PrefixIncrement, :PrefixDecrement,
           :AddressOperator, :IndirectionOperator,
           :One2OneMappedUnaryOperator]
    @eval struct $op <: UnaryOperator end
end

const C2JULIA_UNARY_OPERATOR_MAP = Dict("+" => One2OneMappedUnaryOperator(),
                                        "-" => One2OneMappedUnaryOperator(),
                                        "!" => One2OneMappedUnaryOperator(),
                                        "~" => One2OneMappedUnaryOperator(),
                                        "sizeof" => One2OneMappedUnaryOperator(),
                                        "++" => PrefixIncrement(),
                                        "--" => PrefixDecrement(),
                                        "&" => AddressOperator(),
                                        "*" => IndirectionOperator(),
                                        )

function translate(::PrefixIncrement, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    if typeof(toks[1]) == Punctuation
        return Expr(:macrocall, Symbol("@+"), nothing, Symbol(toks[2].text))
    else
        return Expr(:(+=), Symbol(toks[1].text), 1)
    end
end

function translate(::PrefixDecrement, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    if typeof(toks[1]) == Punctuation
        return Expr(:macrocall, Symbol("@-"), nothing, Symbol(toks[2].text))
    else
        return Expr(:(-=), Symbol(toks[1].text), 1)
    end
end

function translate(::One2OneMappedUnaryOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    operand = children(cursor) |> first |> translate
    Expr(:call, Symbol(op), operand)
end

translate(cursor::CLCStyleCastExpr) = Expr(:call, Symbol(clang2julia(type(cursor))), translate(children(cursor)[]))

function translate(cursor::CLUnaryOperator)
    op = ""
    toks = tokenize(cursor)
    for tok in toks
        typeof(tok) == Punctuation && (op = tok.text)
    end
    return translate(C2JULIA_UNARY_OPERATOR_MAP[op], cursor, op, toks)
end

translate(::AddressOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList) = "translate rule for unary operator & has not implemented yet"
translate(::IndirectionOperator, cursor::CLUnaryOperator) = "translate rule for unary operator * has not implemented yet"



C2JULIA_BINARY_OPERATOR_MAP = Dict()

function translate(cursor::CLBinaryOperator)
    child_cursors = children(cursor)
    toks = tokenize(cursor)
    op = toks[2].text
    if haskey(C2JULIA_BINARY_OPERATOR_MAP, op)
        op_sym = C2JULIA_BINARY_OPERATOR_MAP[op]
    else
        op_sym = Symbol(op)
    end
    lhs = first(child_cursors) |> translate
    rhs = last(child_cursors) |> translate
    return Expr(:call, op_sym, lhs, rhs)
end
