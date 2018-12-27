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

function translate(::AddressOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    return children(cursor) |> first |> translate
end

function translate(::IndirectionOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    return children(cursor) |> first |> translate
end

function translate(::One2OneMappedUnaryOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    operand = children(cursor) |> first |> translate
    Expr(:call, Symbol(op), operand)
end

function translate(cursor::CLUnaryOperator)
    op = ""
    toks = tokenize(cursor)
    # hmm, this is just a workaround, pending https://reviews.llvm.org/D10833
    for tok in toks
        txt = tok.text
        if typeof(tok) == Punctuation && haskey(C2JULIA_UNARY_OPERATOR_MAP, txt)
            op = txt
            break
        end
    end
    return translate(C2JULIA_UNARY_OPERATOR_MAP[op], cursor, op, toks)
end

abstract type BinaryOperator end

for op in [:ArithmeticBinaryOperator, :RelationalBinaryOperator,
           :BitwiseBinaryOperator, :LogicalBinaryOperator, :AssignmentBinaryOperator]
    @eval struct $op <: BinaryOperator end
end

const C2JULIA_BINARY_OPERATOR_MAP = Dict("+" => ArithmeticBinaryOperator(),
                                         "-" => ArithmeticBinaryOperator(),
                                         "*" => ArithmeticBinaryOperator(),
                                         "/" => ArithmeticBinaryOperator(),
                                         "%" => ArithmeticBinaryOperator(),
                                         "<" => RelationalBinaryOperator(),
                                         ">" => RelationalBinaryOperator(),
                                         "==" => RelationalBinaryOperator(),
                                         "!=" => RelationalBinaryOperator(),
                                         ">=" => RelationalBinaryOperator(),
                                         "<=" => RelationalBinaryOperator(),
                                         "~" => BitwiseBinaryOperator(),
                                         "^" => BitwiseBinaryOperator(),
                                         "|" => BitwiseBinaryOperator(),
                                         "&" => BitwiseBinaryOperator(),
                                         "<<" => BitwiseBinaryOperator(),
                                         ">>" => BitwiseBinaryOperator(),
                                         "||" => LogicalBinaryOperator(),
                                         "&&" => LogicalBinaryOperator(),
                                         "=" => AssignmentBinaryOperator(),
                                         "+=" => AssignmentBinaryOperator(),
                                         "-=" => AssignmentBinaryOperator(),
                                         "*=" => AssignmentBinaryOperator(),
                                         "/=" => AssignmentBinaryOperator(),
                                         "%=" => AssignmentBinaryOperator(),
                                        )

function translate(cursor::CLBinaryOperator)
    op = ""
    toks = tokenize(cursor)
    # hmm, this is just a workaround, pending https://reviews.llvm.org/D10833
    for tok in toks
        txt = tok.text
        if typeof(tok) == Punctuation && haskey(C2JULIA_BINARY_OPERATOR_MAP, txt)
            op = txt
            if op != "*"  # possible x*
                break
            end
        end
    end
    child_cursors = children(cursor)
    lhs = first(child_cursors) |> translate
    rhs = last(child_cursors) |> translate
    return Expr(:call, Symbol(op), lhs, rhs)
end
