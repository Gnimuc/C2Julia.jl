function translate(cursor::CLDeclStmt)
    child_cursors = children(cursor)
    decl_exprs = MetaExpr[]
    for c in child_cursors
        meta = translate(c)
        if meta.expr isa Symbol
            meta.expr = Expr(:local, meta.expr)
            push!(decl_exprs, meta)
        else
            push!(decl_exprs, meta)
        end
    end
    return decl_exprs
end

function translate(cursor::CLCompoundStmt)
    child_cursors = children(cursor)
    meta = MetaExpr[]
    block = Expr(:block)
    for c in child_cursors
        append!(meta, translate(c))
    end
    for m in meta
        push!(block.args, m.expr)
        # TODO: add space between stmt
    end
    return MetaExpr(block)
end

function translate(cursor::CLIfStmt)
    child_cursors = children(cursor)
    if_expr = Expr(:if)
    for c in child_cursors
        if kind(c) == CXCursor_IfStmt
            elseif_expr = translate(c).expr
            elseif_expr.head = :elseif
            push!(if_expr.args, elseif_expr)
        else
            push!(if_expr.args, translate(c).expr)
        end
    end
    return MetaExpr(if_expr)
end

translate(cursor::CLBreakStmt) = MetaExpr(Expr(:break), cursor)
translate(cursor::CLContinueStmt) = MetaExpr(Expr(:continue), cursor)

function translate(cursor::CLGotoStmt)
    meta = translate(first(children(cursor)))
    meta.expr = Expr(:macrocall, Symbol("@goto"), nothing, meta.expr)
    return meta
end

function translate(cursor::CLLabelStmt)
    ex = Expr(:macrocall, Symbol("@label"), nothing, Symbol(spelling(cursor)))
    return MetaExpr(ex, cursor)
end

translate(cursor::CLLabelRef) = MetaExpr(Symbol(spelling(cursor)), cursor)

function translate(cursor::CLWhileStmt)
    child_cursors = children(cursor)
    condition_meta = translate(first(child_cursors))
    body_meta = translate(last(child_cursors))
    return MetaExpr(Expr(:while, condition_meta.expr, body_meta.expr))
end

function translate(cursor::CLDoStmt)
    child_cursors = children(cursor)
    condition_meta = translate(last(child_cursors))
    body_meta = translate(first(child_cursors))
    push!(body_meta.expr.args, Expr(:call, :||, condition_meta.expr, Expr(:break)))
    return MetaExpr(Expr(:while, true, body_meta.expr))
end

# this is a workaround since libclang somehow striped those null objects in `for( ; ; )`
function get_init_cond_inc_body(cursor::CLForStmt)
    child_cursors = children(cursor)
    toks = tokenize(cursor) |> collect
    left_paren_idx = 2
    right_paren_idx = findfirst(x->x.text==")", toks)
    idxs = findall(x->x.text==";", toks[1:right_paren_idx])
    @show toks
    idx_count = 1
    if idxs[1] == left_paren_idx+1
        init = MetaExpr(Expr(:null), cursor)
    else
        init = translate(child_cursors[idx_count])
        idx_count += 1
    end

    if idxs[2] == idxs[1]+1
        condition = MetaExpr(Expr(:null), cursor)
    else
        condition = translate(child_cursors[idx_count])
        idx_count += 1
    end

    if right_paren_idx == idxs[2]+1
        increment = MetaExpr(Expr(:null), cursor)
    else
        increment = translate(child_cursors[idx_count])
        idx_count += 1
    end

    @assert idx_count == length(child_cursors) "unknown for-loop structure: $cursor"
    body = translate(last(child_cursors))

    return init, condition, increment, body
end

function translate(cursor::CLForStmt)
    init, cond, inc, body = get_init_cond_inc_body(cursor)
    if init isa Vector{MetaExpr}
        init = init[]
    end
    for_expr = Expr(:macrocall, Symbol("@cfor"), nothing, init.expr, cond.expr, inc.expr, body.expr)
    return MetaExpr(for_expr)
end

function translate(cursor::CLReturnStmt)
    child_cursors = children(cursor)
    if isempty(child_cursors)
        return MetaExpr(Expr(:return, nothing), cursor)
    else
        meta = translate(first(child_cursors))
        meta.expr = Expr(:return, meta.expr)
        return meta
    end
end

translate(cursor::CLNullStmt) = MetaExpr(Expr(:null), cursor)

function translate(cursor::CLCaseStmt)
    child_cursors = children(cursor)
    case_meta = translate(first(child_cursors))
    case_meta.expr = Expr(:macrocall, Symbol("@case"), nothing, case_meta.expr)
    exprs = MetaExpr[case_meta]
    append!(exprs, translate(last(child_cursors)))
    return exprs
end

function translate(cursor::CLDefaultStmt)
    child_cursors = children(cursor)
    default_meta = MetaExpr(Expr(:macrocall, Symbol("@default"), nothing))
    exprs = MetaExpr[default_meta]
    append!(exprs, translate(last(child_cursors)))
    return exprs
end

function translate(cursor::CLSwitchStmt)
    child_cursors = children(cursor)
    constexpr_meta = translate(first(child_cursors))
    body_meta = translate(last(child_cursors))
    ex = Expr(:macrocall, Symbol("@switch"), nothing, constexpr_meta.expr, body_meta.expr)
    return MetaExpr(ex)
end
