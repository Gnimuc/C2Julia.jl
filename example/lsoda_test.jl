function check_opt(ctx, opt)
    mxstp0 = 500
    mord = [0, 12, 5]
    if ctx.state == 0
        ctx.state = 1
    end
    if ctx.state == 1
        opt.h0 = 0.0
        opt.mxordn = mord[2]
        opt.mxords = mord[3]
    end
    if ctx.neq <= 0
        @warn "ERROR(\"[lsoda] neq = %d is less than 1\\n\",ctx->neq)"
        return 0
    end
    rtol = opt.rtol - 1
    atol = opt.atol - 1
    if ctx.state == 1 || ctx.state == 3
        local i
        @cfor i = 1 i <= ctx.neq i += 1 begin
                rtoli = rtol[i + 1]
                atoli = atol[i + 1]
                if rtoli < 0.0
                    @warn "ERROR(\"[lsoda] rtol = %g is less than 0.\\n\",rtoli)"
                end
                if atoli < 0.0
                    @warn "ERROR(\"[lsoda] atol = %g is less than 0.\\n\",atoli)"
                    return 0
                end
            end
    end
    if opt.itask == 0
        opt.itask = 1
    end
    if opt.itask < 1 || opt.itask > 5
        fprintf(__stderrp, "[lsoda] illegal itask = %d\n", opt.itask)
        return 0
    end
    if opt.ixpr < 0 || opt.ixpr > 1
        fprintf(__stderrp, "[lsoda] ixpr = %d is illegal\n", opt.ixpr)
        return 0
    end
    if opt.mxstep < 0
        fprintf(__stderrp, "[lsoda] mxstep < 0\n")
        return 0
    end
    if opt.mxstep == 0
        opt.mxstep = mxstp0
    end
    if opt.mxhnil < 0
        fprintf(__stderrp, "[lsoda] mxhnil < 0\n")
        return 0
    end
    if ctx.state == 1
        if opt.mxordn < 0
            fprintf(__stderrp, "[lsoda] mxordn = %d is less than 0\n", opt.mxordn)
            return 0
        end
        if opt.mxordn == 0
            opt.mxordn = 100
        end
        opt.mxordn = @warn("min(opt->mxordn,mord[1])")
        if opt.mxords < 0
            fprintf(__stderrp, "[lsoda] mxords = %d is less than 0\n", opt.mxords)
            return 0
        end
        if opt.mxords == 0
            opt.mxords = 100
        end
        opt.mxords = @warn("min(opt->mxords,mord[2])")
    end
    if opt.hmax < 0.0
        fprintf(__stderrp, "[lsoda] hmax < 0.\n")
        return 0
    end
    opt.hmxi = 0.0
    if opt.hmax > 0
        opt.hmxi = 1.0 / opt.hmax
    end
    if opt.hmin < 0.0
        fprintf(__stderrp, "[lsoda] hmin < 0.\n")
        return 0
    end
    return 1
end
