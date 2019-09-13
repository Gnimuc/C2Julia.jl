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
        begin
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
function alloc_mem(ctx)
    nyh = ctx.neq
    lenyh = 1 + @warn("max(ctx->opt->mxordn,ctx->opt->mxords)")
    offset = 0
    local i
    yhoff = offset
    offset += (1 + lenyh) * nothing
    yh0off = offset
    begin
        offset += (1 + nyh) * nothing
    end
    wmoff = offset
    wm0off = offset
    offset += (1 + nyh) * nothing
    begin
        offset += (1 + nyh) * nothing
    end
    ewtoff = offset
    offset += (1 + nyh) * nothing
    savfoff = offset
    offset += (1 + nyh) * nothing
    acoroff = offset
    offset += (1 + nyh) * nothing
    ipvtoff = offset
    offset += (1 + nyh) * nothing
    @warn("_C(memory)") = malloc(offset)
    @warn("_C(yh)") = Ptr{Ptr{Cdouble}}(Cstring(@warn("_C(memory)")) + yhoff)
    @warn("_C(wm)") = Ptr{Ptr{Cdouble}}(Cstring(@warn("_C(memory)")) + wmoff)
    @warn("_C(ewt)") = Ptr{Cdouble}(Cstring(@warn("_C(memory)")) + ewtoff)
    @warn("_C(savf)") = Ptr{Cdouble}(Cstring(@warn("_C(memory)")) + savfoff)
    @warn("_C(acor)") = Ptr{Cdouble}(Cstring(@warn("_C(memory)")) + acoroff)
    @warn("_C(ipvt)") = Ptr{Cint}(Cstring(@warn("_C(memory)")) + ipvtoff)
    begin
        (@warn("_C(yh)"))[i + 1] = Ptr{Cdouble}((Cstring(@warn("_C(memory)")) + yh0off) + (i * (1 + nyh)) * nothing)
    end
    begin
        (@warn("_C(wm)"))[i + 1] = Ptr{Cdouble}((Cstring(@warn("_C(memory)")) + wm0off) + (i * (1 + nyh)) * nothing)
    end
    return @warn("_C(memory)") != @warn("NULL")
end
function lsoda_prepare(ctx, opt)
    ctx.common = calloc(1, lsoda_common_t())
    ctx.opt = opt
    if !(check_opt(ctx, opt))
        return 0
    end
    return alloc_mem(ctx)
end
function lsoda_reset(ctx)
    offset = lsoda_common_t() + nothing
    __builtin___memset_chk((@warn("memset"), offset), memset, (lsoda_common_t(), offset), __builtin_object_size((@warn("memset"), offset), memset))
end
function lsoda_free(ctx)
    free((ctx.common).memory)
    if ctx.error
        fprintf(__stderrp, "unhandled error message: %s\n", ctx.error)
        free(ctx.error)
    end
    free(ctx.common)
end
function lsoda(ctx, y, t, tout)
    local kflag
    local jstart
    common = lsoda_common_t(ctx, :common)
    opt = lsoda_opt_t(ctx, :opt)
    y -= 1
    i = 0
    local ihit
    neq = ctx.neq
    local big
    local h0
    local hmx
    local rh
    local tcrit
    local tdist
    local tnext
    local tol
    local tolsf
    local tp
    local size
    local sum
    local w0
    if common == @warn("NULL")
        begin
            @warn "hardfailure(\"[lsoda] illegal common block did you call lsoda_prepare?%s\\n\",\"\")"
            (ctx.state, ?_?(hardfailure))
            return ctx.state
        end
        nothing
    end
    if ctx.state == 1 || ctx.state == 3
        h0 = opt.h0
        if ctx.state == 1
            if (tout - t[]) * h0 < 0.0
                begin
                    @warn "hardfailure(\"[lsoda] tout = %g behind t = %g. integration direction is given by %g\\n\",tout,*t,h0)"
                    (ctx.state, hardfailure[])
                    return ctx.state
                end
                nothing
            end
        end
    end
    itask = opt.itask
    if ctx.state == 3
        jstart = -1
    end
    rtol = opt.rtol - 1
    atol = opt.atol - 1
    if ctx.state == 1
        @warn("_C(meth)") = 1
        @warn("_C(tn)") = t[]
        @warn("_C(tsw)") = t[]
        if itask == 4 || itask == 5
            tcrit = opt.tcrit
            if (tcrit - tout) * (tout - t[]) < 0.0
                begin
                    @warn "hardfailure(\"[lsoda] itask = 4 or 5 and tcrit behind tout%s\\n\",\"\")"
                    (ctx.state, ?_?(hardfailure))
                    return ctx.state
                end
                nothing
            end
            if h0 != 0.0 && ((t[] + h0) - tcrit) * h0 > 0.0
                h0 = tcrit - t[]
            end
        end
        jstart = 0
        @warn("_C(nq)") = 1
        (t[], y + 1, (@warn("_C(yh)"))[3] + 1, ctx.data)
        @warn("_C(nfe)") = 1
        ((@warn("_C(yh)"))[2])[i + 1] = y[i + 1]
        begin
            local i
            ((@warn("ewset(y)"))[i + 1], ((rtol[i + 1], fabs((@warn("ewset(y)"))[i + 1])), atol[i + 1]))
            begin
                ((@warn("ewset(y)"))[i + 1], (ewset, (@warn("ewset(y)"))[i + 1]))
            end
        end
        nothing
        begin
            if (@warn("_C(ewt)"))[i + 1] <= 0.0
                begin
                    @warn "hardfailure(\"[lsoda] ewt[%d] = %g <= 0.\\n\",i,_C(ewt)[i])"
                    (ctx.state, ?_?(hardfailure))
                    return ctx.state
                end
                nothing
            end
        end
        if h0 == 0.0
            tdist = fabs(tout - t[])
            w0 = fmax(fabs(t[]), fabs(tout))
            if tdist < (2.0ETA) * w0
                begin
                    @warn "hardfailure(\"[lsoda] tout too close to t to start integration%s\\n\",\"\")"
                    (ctx.state, ?_?(hardfailure))
                    return ctx.state
                end
                nothing
            end
            tol = 0.0
            nothing
            tol = fmax(tol, rtol[i + 1])
            if tol <= 0.0
                begin
                    atoli = atol[i + 1]
                    ayi = fabs(y[i + 1])
                    if ayi != 0.0
                        tol = fmax(tol, atoli / ayi)
                    end
                end
            end
            tol = fmax(tol, 100.0ETA)
            tol = fmin(tol, 0.001)
            sum = vmnorm(neq, (@warn("_C(yh)"))[3], @warn("_C(ewt)"))
            sum = 1.0 / ((tol * w0) * w0) + (tol * sum) * sum
            h0 = 1.0 / sqrt(sum)
            h0 = fmin(h0, tdist)
            h0 = h0 * if tout - t[] >= 0.0
                        1.0
                    else
                        -1.0
                    end
        end
        rh = fabs(h0) * opt.hmxi
        if rh > 1.0
            h0 /= rh
        end
        @warn("_C(h)") = h0
        ((@warn("_C(yh)"))[3])[i + 1] *= h0
    end
    if ctx.state == 2 || ctx.state == 3
        jstart = 1
        @warn("_C(nslast)") = @warn("_C(nst)")
        @switch itask begin
                @case 1
                if (@warn("_C(tn)") - tout) * @warn("_C(h)") >= 0.0
                    begin
                        iflag = intdy(ctx, tout, intdyreturn, y)
                        if (iflag, intdyreturn)
                            @warn "intdyreturn()"
                            (y[i + 1], ((@warn("intdyreturn()"))[intdyreturn + 1])[i + 1])
                            (?_?(t), @warn("intdyreturn()"))
                        end
                        (?_?(t), tout)
                        (ctx.state, intdyreturn)
                        return ctx.state
                    end
                    nothing
                end
                break
                @case 2
                break
                @case 3
                tp = @warn("_C(tn)") - @warn("_C(hu)") * (1.0 + 100.0ETA)
                if (tp - tout) * @warn("_C(h)") > 0.0
                    begin
                        @warn "hardfailure(\"[lsoda] itask = %d and tout behind tcur - _C(hu)\\n\",itask)"
                        (ctx.state, ?_?(hardfailure))
                        return ctx.state
                    end
                    nothing
                end
                if (@warn("_C(tn)") - tout) * @warn("_C(h)") < 0.0
                    break
                end
                begin
                    local i
                    neq = ctx.neq
                    (y[i + 1], ((@warn("successreturn()"))[successreturn + 1])[i + 1])
                    (?_?(t), @warn("successreturn()"))
                    if ((itask, successreturn), (itask, successreturn))
                        elseif ihit
                            (?_?(t), tcrit)
                        end
                    end
                    (ctx.state, successreturn)
                    return ctx.state
                end
                nothing
                @case 4
                tcrit = opt.tcrit
                if (@warn("_C(tn)") - tcrit) * @warn("_C(h)") > 0.0
                    begin
                        @warn "hardfailure(\"[lsoda] itask = 4 or 5 and tcrit behind tcur%s\\n\",\"\")"
                        (ctx.state, ?_?(hardfailure))
                        return ctx.state
                    end
                    nothing
                end
                if (tcrit - tout) * @warn("_C(h)") < 0.0
                    begin
                        @warn "hardfailure(\"[lsoda] itask = 4 or 5 and tcrit behind tout%s\\n\",\"\")"
                        (ctx.state, ?_?(hardfailure))
                        return ctx.state
                    end
                    nothing
                end
                if (@warn("_C(tn)") - tout) * @warn("_C(h)") >= 0.0
                    begin
                        iflag = intdy(ctx, tout, intdyreturn, y)
                        if (iflag, intdyreturn)
                            @warn "intdyreturn()"
                            (y[i + 1], ((@warn("intdyreturn()"))[intdyreturn + 1])[i + 1])
                            (?_?(t), @warn("intdyreturn()"))
                        end
                        (?_?(t), tout)
                        (ctx.state, intdyreturn)
                        return ctx.state
                    end
                    nothing
                end
                @case 5
                if itask == 5
                    tcrit = opt.tcrit
                    if (@warn("_C(tn)") - tcrit) * @warn("_C(h)") > 0.0
                        begin
                            @warn "hardfailure(\"[lsoda] itask = 4 or 5 and tcrit behind tcur%s\\n\",\"\")"
                            (ctx.state, ?_?(hardfailure))
                            return ctx.state
                        end
                        nothing
                    end
                end
                hmx = fabs(@warn("_C(tn)")) + fabs(@warn("_C(h)"))
                ihit = fabs(@warn("_C(tn)") - tcrit) <= (100.0ETA) * hmx
                if ihit
                    t[] = tcrit
                    begin
                        local i
                        neq = ctx.neq
                        (y[i + 1], ((@warn("successreturn()"))[successreturn + 1])[i + 1])
                        (?_?(t), @warn("successreturn()"))
                        if ((itask, successreturn), (itask, successreturn))
                            elseif ihit
                                (?_?(t), tcrit)
                            end
                        end
                        (ctx.state, successreturn)
                        return ctx.state
                    end
                    nothing
                end
                tnext = @warn("_C(tn)") + @warn("_C(h)") * (1.0 + 4.0ETA)
                if (tnext - tcrit) * @warn("_C(h)") <= 0.0
                    break
                end
                @warn("_C(h)") = (tcrit - @warn("_C(tn)")) * (1.0 - 4.0ETA)
                if ctx.state == 2
                    jstart = -2
                end
                break
            end
    end
    while 1
        if ctx.state != 1 || @warn("_C(nst)") != 0
            if @warn("_C(nst)") - @warn("_C(nslast)") >= opt.mxstep
                begin
                    local i
                    neq = ctx.neq
                    @warn "softfailure(-1,\"[lsoda] %d steps taken before reaching tout\\n\",opt->mxstep)"
                    (y[i + 1], ((@warn("softfailure(-1,\"[lsoda] %d steps taken before reaching tout\\n\",opt->mxstep)"))[softfailure + 1])[i + 1])
                    (-t, @warn("softfailure(-1,\"[lsoda] %d steps taken before reaching tout\\n\",opt->mxstep)"))
                    (ctx.state, ?_?(softfailure))
                    return ctx.state
                end
                nothing
            end
            begin
                local i
                ((@warn("ewset(_C(yh)[1])"))[i + 1], ((rtol[i + 1], fabs((@warn("ewset(_C(yh)[1])"))[i + 1])), atol[i + 1]))
                begin
                    ((@warn("ewset(_C(yh)[1])"))[i + 1], (ewset, (@warn("ewset(_C(yh)[1])"))[i + 1]))
                end
            end
            nothing
            begin
                if (@warn("_C(ewt)"))[i + 1] <= 0.0
                    begin
                        local i
                        neq = ctx.neq
                        @warn "softfailure(-6,\"[lsoda] ewt[%d] = %g <= 0.\\n\",i,_C(ewt)[i])"
                        (y[i + 1], ((@warn("softfailure(-6,\"[lsoda] ewt[%d] = %g <= 0.\\n\",i,_C(ewt)[i])"))[softfailure + 1])[i + 1])
                        (-t, @warn("softfailure(-6,\"[lsoda] ewt[%d] = %g <= 0.\\n\",i,_C(ewt)[i])"))
                        (ctx.state, ?_?(softfailure))
                        return ctx.state
                    end
                    nothing
                end
            end
        end
        tolsf = ETA * vmnorm(neq, (@warn("_C(yh)"))[2], @warn("_C(ewt)"))
        if tolsf > 0.01
            tolsf = tolsf * 200.0
            if @warn("_C(nst)") == 0
                begin
                    @warn "hardfailure(\"lsoda -- at start of problem, too much accuracy\\n\"\" requested for precision of machine,\\n\"\" suggested scaling factor = %g\\n\",tolsf)"
                    (ctx.state, ?_?(hardfailure))
                    return ctx.state
                end
                nothing
            end
            begin
                local i
                neq = ctx.neq
                @warn "softfailure(-2,\"lsoda -- at t = %g, too much accuracy requested\\n\"\"         for precision of machine, suggested\\n\"\"         scaling factor = %g\\n\",*t,tolsf)"
                (y[i + 1], ((@warn("softfailure(-2,\"lsoda -- at t = %g, too much accuracy requested\\n\"\"         for precision of machine, suggested\\n\"\"         scaling factor = %g\\n\",*t,tolsf)"))[softfailure + 1])[i + 1])
                (-t, @warn("softfailure(-2,\"lsoda -- at t = %g, too much accuracy requested\\n\"\"         for precision of machine, suggested\\n\"\"         scaling factor = %g\\n\",*t,tolsf)"))
                (ctx.state, ?_?(softfailure))
                return ctx.state
            end
            nothing
        end
        if @warn("_C(tn)") + @warn("_C(h)") == @warn("_C(tn)")
            @warn("_C(nhnil)") += 1
            if @warn("_C(nhnil)") <= opt.mxhnil
                fprintf(__stderrp, "lsoda -- warning..internal t = %g and _C(h) = %g are\n", @warn("_C(tn)"), @warn("_C(h)"))
                fprintf(__stderrp, "         such that in the machine, t + _C(h) = t on the next step\n")
                fprintf(__stderrp, "         solver will continue anyway.\n")
                if @warn("_C(nhnil)") == opt.mxhnil
                    fprintf(__stderrp, "lsoda -- above warning has been issued %d times,\n", @warn("_C(nhnil)"))
                    fprintf(__stderrp, "         it will not be issued again for this problem\n")
                end
            end
        end
        kflag = stoda(ctx, y, jstart)
        if kflag == 0
            jstart = 1
            if @warn("_C(meth)") != @warn("_C(mused)")
                @warn("_C(tsw)") = @warn("_C(tn)")
                jstart = -1
                if opt.ixpr
                    if @warn("_C(meth)") == 2
                        fprintf(__stderrp, "[lsoda] a switch to the stiff method has occurred ")
                    end
                    if @warn("_C(meth)") == 1
                        fprintf(__stderrp, "[lsoda] a switch to the nonstiff method has occurred")
                    end
                    fprintf(__stderrp, "at t = %g, tentative step size _C(h) = %g, step _C(nst) = %d\n", @warn("_C(tn)"), @warn("_C(h)"), @warn("_C(nst)"))
                end
            end
            if itask == 1
                if (@warn("_C(tn)") - tout) * @warn("_C(h)") < 0.0
                    continue
                end
                begin
                    iflag = intdy(ctx, tout, intdyreturn, y)
                    if (iflag, intdyreturn)
                        @warn "intdyreturn()"
                        (y[i + 1], ((@warn("intdyreturn()"))[intdyreturn + 1])[i + 1])
                        (?_?(t), @warn("intdyreturn()"))
                    end
                    (?_?(t), tout)
                    (ctx.state, intdyreturn)
                    return ctx.state
                end
                nothing
            end
            if itask == 2
                begin
                    local i
                    neq = ctx.neq
                    (y[i + 1], ((@warn("successreturn()"))[successreturn + 1])[i + 1])
                    (?_?(t), @warn("successreturn()"))
                    if ((itask, successreturn), (itask, successreturn))
                        elseif ihit
                            (?_?(t), tcrit)
                        end
                    end
                    (ctx.state, successreturn)
                    return ctx.state
                end
                nothing
            end
            if itask == 3
                if (@warn("_C(tn)") - tout) * @warn("_C(h)") >= 0.0
                    begin
                        local i
                        neq = ctx.neq
                        (y[i + 1], ((@warn("successreturn()"))[successreturn + 1])[i + 1])
                        (?_?(t), @warn("successreturn()"))
                        if ((itask, successreturn), (itask, successreturn))
                            elseif ihit
                                (?_?(t), tcrit)
                            end
                        end
                        (ctx.state, successreturn)
                        return ctx.state
                    end
                    nothing
                end
                continue
            end
            if itask == 4
                tcrit = opt.tcrit
                if (@warn("_C(tn)") - tout) * @warn("_C(h)") >= 0.0
                    begin
                        iflag = intdy(ctx, tout, intdyreturn, y)
                        if (iflag, intdyreturn)
                            @warn "intdyreturn()"
                            (y[i + 1], ((@warn("intdyreturn()"))[intdyreturn + 1])[i + 1])
                            (?_?(t), @warn("intdyreturn()"))
                        end
                        (?_?(t), tout)
                        (ctx.state, intdyreturn)
                        return ctx.state
                    end
                    nothing
                else
                    hmx = fabs(@warn("_C(tn)")) + fabs(@warn("_C(h)"))
                    ihit = fabs(@warn("_C(tn)") - tcrit) <= (100.0ETA) * hmx
                    if ihit
                        begin
                            local i
                            neq = ctx.neq
                            (y[i + 1], ((@warn("successreturn()"))[successreturn + 1])[i + 1])
                            (?_?(t), @warn("successreturn()"))
                            if ((itask, successreturn), (itask, successreturn))
                                elseif ihit
                                    (?_?(t), tcrit)
                                end
                            end
                            (ctx.state, successreturn)
                            return ctx.state
                        end
                        nothing
                    end
                    tnext = @warn("_C(tn)") + @warn("_C(h)") * (1.0 + 4.0ETA)
                    if (tnext - tcrit) * @warn("_C(h)") <= 0.0
                        continue
                    end
                    @warn("_C(h)") = (tcrit - @warn("_C(tn)")) * (1.0 - 4.0ETA)
                    jstart = -2
                    continue
                end
            end
            if itask == 5
                tcrit = opt.tcrit
                hmx = fabs(@warn("_C(tn)")) + fabs(@warn("_C(h)"))
                ihit = fabs(@warn("_C(tn)") - tcrit) <= (100.0ETA) * hmx
                begin
                    local i
                    neq = ctx.neq
                    (y[i + 1], ((@warn("successreturn()"))[successreturn + 1])[i + 1])
                    (?_?(t), @warn("successreturn()"))
                    if ((itask, successreturn), (itask, successreturn))
                        elseif ihit
                            (?_?(t), tcrit)
                        end
                    end
                    (ctx.state, successreturn)
                    return ctx.state
                end
                nothing
            end
        end
        if kflag == -1 || kflag == -2
            big = 0.0
            @warn("_C(imxer)") = 1
            begin
                size = fabs((@warn("_C(acor)"))[i + 1]) * (@warn("_C(ewt)"))[i + 1]
                if big < size
                    big = size
                    @warn("_C(imxer)") = i
                end
            end
            if kflag == -1
                begin
                    local i
                    neq = ctx.neq
                    @warn "softfailure(-4,\"lsoda -- at t = %g and step size _C(h) = %g, the\\n\",\"         error test failed repeatedly or\\n\"\"         with fabs(_C(h)) = hmin\\n\",_C(tn),_C(h))"
                    (y[i + 1], ((@warn("softfailure(-4,\"lsoda -- at t = %g and step size _C(h) = %g, the\\n\",\"         error test failed repeatedly or\\n\"\"         with fabs(_C(h)) = hmin\\n\",_C(tn),_C(h))"))[softfailure + 1])[i + 1])
                    (-t, @warn("softfailure(-4,\"lsoda -- at t = %g and step size _C(h) = %g, the\\n\",\"         error test failed repeatedly or\\n\"\"         with fabs(_C(h)) = hmin\\n\",_C(tn),_C(h))"))
                    (ctx.state, ?_?(softfailure))
                    return ctx.state
                end
                nothing
            end
            if kflag == -2
                begin
                    local i
                    neq = ctx.neq
                    @warn "softfailure(-5,\"lsoda -- at t = %g and step size _C(h) = %g, the\\n\"\"         corrector convergence failed repeatedly or\\n\"\"         with fabs(_C(h)) = hmin\\n\",_C(tn),_C(h))"
                    (y[i + 1], ((@warn("softfailure(-5,\"lsoda -- at t = %g and step size _C(h) = %g, the\\n\"\"         corrector convergence failed repeatedly or\\n\"\"         with fabs(_C(h)) = hmin\\n\",_C(tn),_C(h))"))[softfailure + 1])[i + 1])
                    (-t, @warn("softfailure(-5,\"lsoda -- at t = %g and step size _C(h) = %g, the\\n\"\"         corrector convergence failed repeatedly or\\n\"\"         with fabs(_C(h)) = hmin\\n\",_C(tn),_C(h))"))
                    (ctx.state, ?_?(softfailure))
                    return ctx.state
                end
                nothing
            end
        end
    end
end
function lsoda_create_ctx()
    mem = lsoda_context_t(malloc, lsoda_context_t())
    return mem
end
function lsoda_create_opt()
    mem = lsoda_opt_t(malloc, lsoda_opt_t())
    return mem
end
function lsoda_free_opt(opt)
    free(opt.atol)
    free(opt.rtol)
    free(opt)
end
