module CSyntax

export @+, @-

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

end  # module
