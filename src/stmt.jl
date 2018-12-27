function translate(cursor::CLDeclStmt)::Vector{Expr}
    child_cursors = children(cursor)
    decl_exprs = Expr[]
    for c in child_cursors
        subexpr = translate(c)
        if subexpr isa Symbol
            push!(decl_exprs, Expr(:local, translate(c)))
        else
            push!(decl_exprs, subexpr)
        end
    end
    return decl_exprs
end

function translate(cursor::CLCompoundStmt)
    child_cursors = children(cursor)
    block = Expr(:block)
    for c in child_cursors
        expr = translate(c)
        if expr isa Vector
            push!(block.args, expr...)
        else
            push!(block.args, expr)
        end
    end
    return block
end

function translate(cursor::CLIfStmt)
    child_cursors = children(cursor)
    if_expr = Expr(:if)
    for c in child_cursors
        if kind(c) == CXCursor_IfStmt
            elseif_expr = translate(c)
            elseif_expr.head = :elseif
            push!(if_expr.args, elseif_expr)
        else
            push!(if_expr.args, translate(c))
        end
    end
    return if_expr
end

translate(cursor::CLBreakStmt) = Expr(:break)
translate(cursor::CLContinueStmt) = Expr(:continue)
translate(cursor::CLGotoStmt) = Expr(:macrocall, Symbol("@goto"), nothing, translate(children(cursor)[]))
translate(cursor::CLLabelStmt) = Expr(:macrocall, Symbol("@label"), nothing, Symbol(spelling(cursor)))
translate(cursor::CLLabelRef) = Symbol(spelling(cursor))

function translate(cursor::CLWhileStmt)
    child_cursors = children(cursor)
    condition = first(child_cursors) |> translate
    body = last(child_cursors) |> translate
    return Expr(:while, condition, body)
end

function translate(cursor::CLDoStmt)
    child_cursors = children(cursor)
    condition = last(child_cursors) |> translate
    body = first(child_cursors) |> translate
    push!(body.args, Expr(:call, :||, condition, Expr(:break)))
    return Expr(:while, true, body)
end

# this is a workaround since libclang somehow striped those null objects in `for( ; ; )`
function get_init_cond_inc_body(cursor::CLForStmt)
    child_cursors = children(cursor)
    toks = tokenize(cursor) |> collect
    left_paren_idx = 2
    right_paren_idx = findfirst(x->x.text==")", toks)
    idxs = findall(x->x.text==";", toks[1:right_paren_idx])

    idx_count = 1
    if idxs[1] == left_paren_idx+1
        init = Expr(:null)
    else
        init = translate(child_cursors[idx_count])
        idx_count += 1
    end

    if idxs[2] == idxs[1]+1
        condition = Expr(:null)
    else
        condition = translate(child_cursors[idx_count])
        idx_count += 1
    end

    if right_paren_idx == idxs[2]+1
        increment = Expr(:null)
    else
        increment = translate(child_cursors[idx_count])
        idx_count += 1
    end

    @assert idx_count == length(child_cursors) "unknown for-loop structure: $cursor"
    body = last(child_cursors) |> translate

    return init, condition, increment, body
end

function translate(cursor::CLForStmt)
    init, cond, inc, body = get_init_cond_inc_body(cursor)
    let_expr = Expr(:let)
    if init == Expr(:null)
        push!(let_expr.args, Expr(:block))
    else
        push!(let_expr.args, init)
    end
    cond == Expr(:null) && (cond = true;)
    inc != Expr(:null) && push!(body.args, inc)
    push!(let_expr.args, Expr(:while, cond, body))
    return let_expr
end

function translate(cursor::CLReturnStmt)
    child_cursors = children(cursor)
    if isempty(child_cursors)
        Expr(:return, nothing)
    else
        Expr(:return, translate(child_cursors[]))
    end
end

translate(cursor::CLNullStmt) = Expr(:null)

Base.length(::Expr) = 1
Base.iterate(x::Expr) = (x, nothing)

function translate(cursor::CLCaseStmt)::Vector{Expr}
    child_cursors = children(cursor)
    case_expr = Expr(:macrocall, Symbol("@case"), nothing, translate(first(child_cursors)))
    exprs = Expr[case_expr]
    append!(exprs, translate(last(child_cursors)))
    return exprs
end

function translate(cursor::CLDefaultStmt)::Vector{Expr}
    child_cursors = children(cursor)
    default_expr = Expr(:macrocall, Symbol("@default"), nothing)
    exprs = Expr[default_expr]
    append!(exprs, translate(last(child_cursors)))
    return exprs
end

function translate(cursor::CLSwitchStmt)
    child_cursors = children(cursor)
    constexpr = first(child_cursors) |> translate
    body_expr = last(child_cursors) |> translate
    Expr(:macrocall, Symbol("@switch"), nothing, constexpr, body_expr)
end
