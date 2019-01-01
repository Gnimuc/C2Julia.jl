module CSyntax

const malloc = Base.Libc.malloc
const free = Base.Libc.free
export malloc, free

export @+, @++, @-
export @cfor, @switch

# for ref indexing a string
@inline Base.getindex(s::AbstractString) = s[1]

"""
PrefixIncrement Operator
"""
macro +(x)
    @gensym tmp
    quote
        tmp = $(esc(x))
        $(esc(x)) += 1
        tmp
    end
end
@eval const $(Symbol("@++")) = $(Symbol("@+"))

"""
PrefixDecrement Operator
"""
macro -(x)
    @gensym tmp
    quote
        tmp = $(esc(x))
        $(esc(x)) -= 1
        tmp
    end
end
# @-- is invalid

macro cfor(init, cond, inc, body)
    let_expr = Expr(:let, Expr(:block))
    if !(Meta.isexpr(init, :block) && isempty(init.args) || Meta.isexpr(init, :null))
        push!(let_expr.args[1].args, init)
    end
    condition = cond
    if Meta.isexpr(condition, :block) && isempty(condition.args) || Meta.isexpr(init, :null)
        condition = true
    end
    Meta.isexpr(inc, :block) && isempty(inc.args) || Meta.isexpr(init, :null) || push!(body.args, inc)
    push!(let_expr.args, Expr(:while, condition, body))
    return esc(let_expr)
end

# Inspired by dcjones's Switch.jl
macro switch(constexpr, body)
    case2label = Dict{Any,Symbol}()
    flow = Expr(:block)
    end_label = gensym("end")
    default_label = end_label

    for arg in body.args
        if Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@case")
            label = gensym("case")
            case2label[arg.args[3]] = label
            labelexpr = Expr(:symboliclabel, label)
            push!(flow.args, labelexpr)
        elseif Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@default")
            default_label = gensym("default")
            labelexpr = Expr(:symboliclabel, default_label)
            push!(flow.args, labelexpr)
        elseif arg == Expr(:break)
            labelexpr = Expr(:symbolicgoto, end_label)
            push!(flow.args, labelexpr)
        else
            push!(flow.args, arg)
        end
    end
    push!(flow.args, Expr(:symboliclabel, end_label))

    jumptable = Expr(:block)
    for (case, label) in case2label
        condition = Expr(:call, :(==), constexpr, case)
        push!(jumptable.args, Expr(:if, condition, Expr(:symbolicgoto, label)))
    end
    push!(jumptable.args[end].args, Expr(:symbolicgoto, default_label))

    return esc(Expr(:block, jumptable, flow))
end



end  # module
