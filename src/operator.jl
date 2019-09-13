abstract type UnaryOperator end

for op in [:PrefixIncrement, :PrefixDecrement,
           :AddressOperator, :IndirectionOperator,
           :One2OneMappedUnaryOperator, :UnknownOperator]
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
                                        "" => UnknownOperator()
                                        )

function translate(::PrefixIncrement, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    oprand_cursor = children(cursor)[]
    if typeof(first(toks)) == Punctuation
        ex = Expr(:macrocall, Symbol("@++"), nothing, translate(oprand_cursor).expr)
        return MetaExpr(ex, cursor)
    else
        ex = Expr(:(+=), translate(oprand_cursor).expr, 1)
        return MetaExpr(ex, cursor)
    end
end

function translate(::PrefixDecrement, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    oprand_cursor = children(cursor)[]
    if typeof(first(toks)) == Punctuation
        ex = Expr(:macrocall, Symbol("@-"), nothing, translate(oprand_cursor).expr)
        return MetaExpr(ex, cursor)
    else
        ex = Expr(:(-=), translate(oprand_cursor).expr, 1)
        return MetaExpr(ex, cursor)
    end
end

function translate(::AddressOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    meta = translate(first(children(cursor)))
    jltype_sym = clang2julia(type(meta.leafcursor))
    if jltype_sym isa Expr
        # dirty workaround
        jltype_sym = jltype_sym.args[2]
    end
    # if isdefined(Base, jltype_sym) && isbitstype(getfield(Base, jltype_sym))  # TODO: ismutable?
        # immutable isbitstypes needs to be translated to `Ref`s to match C's behavior
        meta.expr = Expr(:call, :Ref, meta.expr)
    # end
    meta.info = "AddressOperator"
    return meta
end

function translate(::IndirectionOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    meta = translate(first(children(cursor)))
    jltype_sym = clang2julia(type(meta.leafcursor))
    if jltype_sym isa Expr
        # dirty workaround
        jltype_sym = jltype_sym.args[2]
    end
    # if isdefined(Base, jltype_sym) && isbitstype(getfield(Base, jltype_sym))  # TODO: ismutable?
        meta.expr = Expr(:ref, meta.expr)
    # end
    return meta
end

function translate(::One2OneMappedUnaryOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    operand_meta = translate(first(children(cursor)))
    return MetaExpr(Expr(:call, Symbol(op), operand_meta.expr), cursor)
end


function translate(::UnknownOperator, cursor::CLUnaryOperator, op::AbstractString, toks::TokenList)
    operand_meta = translate(first(children(cursor)))
    return MetaExpr(Expr(:call, Symbol("?_?"), operand_meta.expr), cursor)
end

function translate(cursor::CLUnaryOperator)
    op = ""
    toks = tokenize(cursor)
    # hmm, this is just a workaround, pending https://reviews.llvm.org/D10833
    GC.@preserve toks begin
        for tok in toks
            txt = tok.text
            if typeof(tok) == Punctuation && haskey(C2JULIA_UNARY_OPERATOR_MAP, txt)
                op = txt
                break
            end
        end
    end
    if op == ""
        @warn "failed to extract unary operator for $cursor !!!"
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
                                         "|=" => AssignmentBinaryOperator(),
                                         "&=" => AssignmentBinaryOperator(),
                                         "<<=" => AssignmentBinaryOperator(),
                                         ">>=" => AssignmentBinaryOperator(),
                                        )

function count_parens(cursor)
    n = 0
    toks = tokenize(cursor)
    for tok in toks
        txt = tok.text
        if txt == "("
            n += 1
        else
            break
        end
    end
    return n
end

function translate(cursor::Union{CLBinaryOperator,CLCompoundAssignOperator})
    # hmm, this is just a workaround, pending https://reviews.llvm.org/D10833
    paren_num = count_parens(cursor)
    child_cursors = children(cursor)
    lhs_cursor = first(child_cursors)
    rhs_cursor = last(child_cursors)
    tmp_cursor = lhs_cursor
    for i = 1:paren_num-1
        xxx = children(tmp_cursor)
        tmp_cursor = first(children(tmp_cursor))
    end
    GC.@preserve tmp_cursor begin
        lhs_last_idx = length(tokenize(tmp_cursor))
        toks = tokenize(cursor)
        op = ""
        for tok in collect(toks)[lhs_last_idx:end]
            txt = tok.text
            if typeof(tok) == Punctuation && haskey(C2JULIA_BINARY_OPERATOR_MAP, txt)
                op = txt
                break
            end
        end
    end
    lhs_meta = translate(lhs_cursor)
    rhs_meta = translate(rhs_cursor)
    return MetaExpr(Expr(:call, Symbol(op), lhs_meta.expr, rhs_meta.expr))
end

function translate(cursor::CLConditionalOperator)
    child_cursors = children(cursor)
    condition_meta = translate(child_cursors[1])
    branch1_meta = translate(child_cursors[2])
    branch2_meta = translate(child_cursors[3])
    ex = Expr(:if, condition_meta.expr, branch1_meta.expr, branch2_meta.expr)
    return MetaExpr(ex)
end
