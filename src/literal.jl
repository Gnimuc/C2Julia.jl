# copied from Clang.jl
# normalize literal with a size suffix
function literally(tok)
    # note: put multi-character first, or it will break out too soon for those!
    literalsuffixes = ["ULL", "Ull", "uLL", "ull", "LLU", "LLu", "llU", "llu",
                       "LL", "ll", "UL", "Ul", "uL", "ul", "LU", "Lu", "lU", "lu",
                       "U", "u", "L", "l", "F", "f"]

    function literal_totype(literal, txt)
      literal = lowercase(literal)

      # Floats following http://en.cppreference.com/w/cpp/language/floating_literal
      float64 = occursin(".", txt) && occursin("l", literal)
      float32 = occursin("f", literal)

      if float64 || float32
        float64 && return "Float64"
        float32 && return "Float32"
      end

      # Integers following http://en.cppreference.com/w/cpp/language/integer_literal
      unsigned = occursin("u", literal)
      nbits = count(x -> x == 'l', literal) == 2 ? 64 : 32
      return "$(unsigned ? "U" : "")Int$nbits"
    end

    token_kind = kind(tok)
    txt = tok.text |> strip
    if token_kind == CXToken_Identifier || token_kind == CXToken_Punctuation
        # pass
    elseif token_kind == CXToken_Literal
        for sfx in literalsuffixes
            if endswith(txt, sfx)
                type = literal_totype(sfx, txt)
                txt = txt[1:end-length(sfx)]
                txt = "$(type)($txt)"
                break
            end
        end
    end
    return txt
end


function translate(cursor::CLIntegerLiteral)
    tok = tokenize(cursor)[1]
    if tok.text != "0" &&
       startswith(tok.text, "0") &&
       !startswith(lowercase(tok.text), "0x")
        ex = Meta.parse("0o"*literally(tok))
        return MetaExpr(ex, cursor)
    else
        ex = Meta.parse(literally(tok))
        return MetaExpr(ex, cursor)
    end
end

function translate(cursor::CLFloatingLiteral)
    ex = tokenize(cursor)[1] |> literally |> Meta.parse
    return MetaExpr(ex, cursor)
end

function translate(cursor::CLCharacterLiteral)
    ex = tokenize(cursor)[1].text |> Meta.parse
    return MetaExpr(ex, cursor)
end

function translate(cursor::CLStringLiteral)
    ex = tokenize(cursor)[1].text |> Meta.parse
    return MetaExpr(ex, cursor)
end
