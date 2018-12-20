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
# handy macro for ++x
macro plusplus(x)
    @gensym tmp
    quote
        tmp = $(esc(x))
        $(esc(x)) += 1
        tmp
    end
end

# handy macro for --x
macro minusminus(x)
    @gensym tmp
    quote
        tmp = $(esc(x))
        $(esc(x)) -= 1
        tmp
    end
end

function translate(::PrefixIncrement, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    if typeof(toks[1]) == Punctuation
        return Expr(:macrocall, Symbol("@plusplus"), nothing, Symbol(toks[2].text))
    else
        return Expr(:(+=), Symbol(toks[1].text), 1)
    end
end

function translate(::PrefixDecrement, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    if typeof(toks[1]) == Punctuation
        return Expr(:macrocall, Symbol("@minusminus"), nothing, Symbol(toks[2].text))
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
