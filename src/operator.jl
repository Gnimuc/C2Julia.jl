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
        ex = Expr(:macrocall, Symbol("@+"), nothing, Symbol(toks[2].text))
        return MetaExpr(ex, cursor)
    else
        ex = Expr(:(+=), Symbol(toks[1].text), 1)
        return MetaExpr(ex, cursor)
    end
end

function translate(::PrefixDecrement, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    if typeof(toks[1]) == Punctuation
        ex = Expr(:macrocall, Symbol("@-"), nothing, Symbol(toks[2].text))
        return MetaExpr(ex, cursor)
    else
        ex = Expr(:(-=), Symbol(toks[1].text), 1)
        return MetaExpr(ex, cursor)
    end
end

function translate(::AddressOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    meta = translate(first(children(cursor)))
    jltype_sym = clang2julia(type(meta.leafcursor))
    if isdefined(Base, jltype_sym) && isbitstype(getfield(Base, jltype_sym))  # TODO: ismutable?
        # immutable isbitstypes needs to be translated to `Ref`s to match C's behavior
        meta.expr = Expr(:call, :Ref, meta.expr)
    end
    return meta
end

function translate(::IndirectionOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    meta = translate(first(children(cursor)))
    jltype_sym = clang2julia(type(meta.leafcursor))
    if jltype_sym isa Expr
        # dirty workaround
        jltype_sym = jltype_sym.args[2]
    end
    if isdefined(Base, jltype_sym) && isbitstype(getfield(Base, jltype_sym))  # TODO: ismutable?
        meta.expr = Expr(:ref, meta.expr)
    end
    return meta
end

function translate(::One2OneMappedUnaryOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    operand_meta = translate(first(children(cursor)))
    return MetaExpr(Expr(:call, Symbol(op), operand_meta.expr), cursor)
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
    lhs_meta = translate(first(child_cursors))
    rhs_meta = translate(last(child_cursors))
    return MetaExpr(Expr(:call, Symbol(op), lhs_meta.expr, rhs_meta.expr))
end
