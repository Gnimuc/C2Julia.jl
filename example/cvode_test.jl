function CVAllocVectors(cv_mem, neq, maxord, machEnv)
    local i
    local j
    @warn("ewt") = N_VNew(neq, machEnv)
    if @warn("ewt") == @warn("NULL")
        return FALSE
    end
    @warn("acor") = N_VNew(neq, machEnv)
    if @warn("acor") == @warn("NULL")
        N_VFree(@warn("ewt"))
        return FALSE
    end
    @warn("tempv") = N_VNew(neq, machEnv)
    if @warn("tempv") == @warn("NULL")
        N_VFree(@warn("ewt"))
        N_VFree(@warn("acor"))
        return FALSE
    end
    @warn("ftemp") = N_VNew(neq, machEnv)
    if @warn("ftemp") == @warn("NULL")
        N_VFree(@warn("tempv"))
        N_VFree(@warn("ewt"))
        N_VFree(@warn("acor"))
        return FALSE
    end
    begin
        (@warn("zn"))[j + 1] = N_VNew(neq, machEnv)
        if (@warn("zn"))[j + 1] == @warn("NULL")
            N_VFree(@warn("ewt"))
            N_VFree(@warn("acor"))
            N_VFree(@warn("tempv"))
            N_VFree(@warn("ftemp"))
            N_VFree((@warn("zn"))[i + 1])
            return FALSE
        end
    end
    @warn("lrw") = (maxord + 5) * neq
    @warn("liw") = 0
    return TRUE
end
function CVFreeVectors(cv_mem, maxord)
    local j
    N_VFree(@warn("ewt"))
    N_VFree(@warn("acor"))
    N_VFree(@warn("tempv"))
    N_VFree(@warn("ftemp"))
    N_VFree((@warn("zn"))[j + 1])
end
function CVEwtSet(cv_mem, rtol, atol, tol_type, ycur, neq)
    @switch tol_type begin
            @case SS
            return CVEwtSetSS(cv_mem, rtol, Ptr{real}(atol), ycur, neq)
            @case SV
            return CVEwtSetSV(cv_mem, rtol, N_Vector(atol), ycur, neq)
        end
end
function CVEwtSetSS(cv_mem, rtol, atol, ycur, neq)
    local rtoli
    local atoli
    rtoli = rtol[]
    atoli = atol[]
    N_VAbs(ycur, @warn("tempv"))
    N_VScale(rtoli, @warn("tempv"), @warn("tempv"))
    N_VAddConst(@warn("tempv"), atoli, @warn("tempv"))
    if N_VMin(@warn("tempv")) <= ZERO
        return FALSE
    end
    N_VInv(@warn("tempv"), @warn("ewt"))
    return TRUE
end
function CVEwtSetSV(cv_mem, rtol, atol, ycur, neq)
    local rtoli
    rtoli = rtol[]
    N_VAbs(ycur, @warn("tempv"))
    N_VLinearSum(rtoli, @warn("tempv"), ONE, atol, @warn("tempv"))
    if N_VMin(@warn("tempv")) <= ZERO
        return FALSE
    end
    N_VInv(@warn("tempv"), @warn("ewt"))
    return TRUE
end
function CVHin(cv_mem, tout)
    local sign
    local count
    local tdiff
    local tdist
    local tround
    local hlb
    local hub
    local hg
    local hgs
    local hnew
    local hrat
    local h0
    local yddnrm
    if (tdiff = tout - @warn("tn")) == ZERO
        return FALSE
    end
    sign = if tdiff > ZERO
            1
        else
            -1
        end
    tdist = @warn("ABS(tdiff)")
    tround = @warn("uround") * @warn("MAX(ABS(tn),ABS(tout))")
    if tdist < TWO * tround
        return FALSE
    end
    hlb = HLB_FACTOR * tround
    hub = CVUpperBoundH0(cv_mem, tdist)
    hg = RSqrt(hlb * hub)
    if hub < hlb
        if sign == -1
            hg = -hg
        end
        @warn("h") = hg
        return TRUE
    end
    count = 0
    begin
        hgs = hg * sign
        yddnrm = CVYddNorm(cv_mem, hgs)
        hnew = if (yddnrm * hub) * hub > TWO
                RSqrt(TWO / yddnrm)
            else
                RSqrt(hg * hub)
            end
        count += 1
        if count >= MAX_ITERS
            break
        end
        hrat = hnew / hg
        if hrat > HALF && hrat < TWO
            break
        end
        if count >= 2 && hrat > TWO
            hnew = hg
            break
        end
        hg = hnew
    end
    h0 = H_BIAS * hnew
    if h0 < hlb
        h0 = hlb
    end
    if h0 > hub
        h0 = hub
    end
    if sign == -1
        h0 = -h0
    end
    @warn("h") = h0
    return TRUE
end
function CVUpperBoundH0(cv_mem, tdist)
    local atoli
    local hub_inv
    local hub
    local vectorAtol
    local temp1
    local temp2
    vectorAtol = @warn("itol") == SV
    if !vectorAtol
        atoli = (Ptr{real}(@warn("abstol")))[]
    end
    temp1 = @warn("tempv")
    temp2 = @warn("acor")
    N_VAbs((@warn("zn"))[1], temp1)
    N_VAbs((@warn("zn"))[2], temp2)
    if vectorAtol
        N_VLinearSum(HUB_FACTOR, temp1, ONE, N_Vector(@warn("abstol")), temp1)
    else
        N_VScale(HUB_FACTOR, temp1, temp1)
        N_VAddConst(temp1, atoli, temp1)
    end
    N_VDiv(temp2, temp1, temp1)
    hub_inv = N_VMaxNorm(temp1)
    hub = HUB_FACTOR * tdist
    if hub * hub_inv > ONE
        hub = ONE / hub_inv
    end
    return hub
end
function CVYddNorm(cv_mem, hg)
    local yddnrm
    N_VLinearSum(hg, (@warn("zn"))[2], ONE, (@warn("zn"))[1], @warn("y"))
    (@warn("N"), @warn("tn") + hg, @warn("y"), @warn("tempv"), @warn("f_data"))
    @warn("nfe") += 1
    N_VLinearSum(ONE, @warn("tempv"), -ONE, (@warn("zn"))[2], @warn("tempv"))
    N_VScale(ONE / hg, @warn("tempv"), @warn("tempv"))
    yddnrm = N_VWrmsNorm(@warn("tempv"), @warn("ewt"))
    return yddnrm
end
function CVStep(cv_mem)
    local saved_t
    local dsm
    local ncf
    local nef
    local nflag
    local kflag
    local passed
    saved_t = @warn("tn")
    ncf = (nef = 0)
    nflag = FIRST_CALL
    if @warn("nst") > 0 && @warn("hprime") != @warn("h")
        CVAdjustParams(cv_mem)
    end
    begin
        CVPredict(cv_mem)
        CVSet(cv_mem)
        nflag = CVnls(cv_mem, nflag)
        kflag = @c(CVHandleNFlag(cv_mem, &nflag, saved_t, &ncf))
        if kflag == ?_?(PREDICT_AGAIN)
            continue
        end
        if kflag != DO_ERROR_TEST
            return kflag
        end
        passed = @c(CVDoErrorTest(cv_mem, &nflag, &kflag, saved_t, &nef, &dsm))
        if !passed && kflag == ?_?(REP_ERR_FAIL)
            return kflag
        end
        if passed
            break
        end
    end
    CVCompleteStep(cv_mem)
    CVPrepareNextStep(cv_mem, dsm)
    return SUCCESS_STEP
end
function CVAdjustParams(cv_mem)
    if @warn("qprime") != @warn("q")
        CVAdjustOrder(cv_mem, @warn("qprime") - @warn("q"))
        @warn("q") = @warn("qprime")
        @warn("L") = @warn("q") + 1
        @warn("qwait") = @warn("L")
    end
    CVRescale(cv_mem)
end
function CVAdjustOrder(cv_mem, deltaq)
    if @warn("q") == 2 && deltaq != 1
        return
    end
    @switch @warn("lmm") begin
            @case ADAMS
            CVAdjustAdams(cv_mem, deltaq)
            break
            @case BDF
            CVAdjustBDF(cv_mem, deltaq)
            break
        end
end
function CVAdjustAdams(cv_mem, deltaq)
    local i
    local j
    local xi
    local hsum
    if deltaq == 1
        N_VConst(ZERO, (@warn("zn"))[@warn("L") + 1])
        return
    end
    (@warn("l"))[i + 1] = ZERO
    (@warn("l"))[2] = ONE
    hsum = ZERO
    begin
        hsum += (@warn("tau"))[j + 1]
        xi = hsum / @warn("hscale")
        (@warn("l"))[i + 1] = (@warn("l"))[i + 1] * xi + (@warn("l"))[(i - 1) + 1]
    end
    (@warn("l"))[(j + 1) + 1] = @warn("q") * ((@warn("l"))[j + 1] / (j + 1))
    N_VLinearSum(-((@warn("l"))[j + 1]), (@warn("zn"))[@warn("q") + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
end
function CVAdjustBDF(cv_mem, deltaq)
    @switch deltaq begin
            @case 1
            CVIncreaseBDF(cv_mem)
            return
            @case -1
            CVDecreaseBDF(cv_mem)
            return
        end
end
function CVIncreaseBDF(cv_mem)
    local alpha0
    local alpha1
    local prod
    local xi
    local xiold
    local hsum
    local A1
    local i
    local j
    (@warn("l"))[i + 1] = ZERO
    (@warn("l"))[3] = (alpha1 = (prod = (xiold = ONE)))
    alpha0 = -ONE
    hsum = @warn("hscale")
    if @warn("q") > 1
        begin
            hsum += (@warn("tau"))[(j + 1) + 1]
            xi = hsum / @warn("hscale")
            prod *= xi
            alpha0 -= ONE / (j + 1)
            alpha1 += ONE / xi
            (@warn("l"))[i + 1] = (@warn("l"))[i + 1] * xiold + (@warn("l"))[(i - 1) + 1]
            xiold = xi
        end
    end
    A1 = (-alpha0 - alpha1) / prod
    N_VScale(A1, (@warn("zn"))[@warn("qmax") + 1], (@warn("zn"))[@warn("L") + 1])
    begin
        N_VLinearSum((@warn("l"))[j + 1], (@warn("zn"))[@warn("L") + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
    end
end
function CVDecreaseBDF(cv_mem)
    local hsum
    local xi
    local i
    local j
    (@warn("l"))[i + 1] = ZERO
    (@warn("l"))[3] = ONE
    hsum = ZERO
    begin
        hsum += (@warn("tau"))[j + 1]
        xi = hsum / @warn("hscale")
        (@warn("l"))[i + 1] = (@warn("l"))[i + 1] * xi + (@warn("l"))[(i - 1) + 1]
    end
    N_VLinearSum(-((@warn("l"))[j + 1]), (@warn("zn"))[@warn("q") + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
end
function CVRescale(cv_mem)
    local j
    local factor
    factor = @warn("eta")
    begin
        N_VScale(factor, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
        factor *= @warn("eta")
    end
    @warn("h") = @warn("hscale") * @warn("eta")
    @warn("hscale") = @warn("h")
end
function CVPredict(cv_mem)
    local j
    local k
    @warn("tn") += @warn("h")
    N_VLinearSum(ONE, (@warn("zn"))[(j - 1) + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[(j - 1) + 1])
end
function CVSet(cv_mem)
    @switch @warn("lmm") begin
            @case ADAMS
            CVSetAdams(cv_mem)
            break
            @case BDF
            CVSetBDF(cv_mem)
            break
        end
    @warn("rl1") = ONE / (@warn("l"))[2]
    @warn("gamma") = @warn("h") * @warn("rl1")
    if @warn("nst") == 0
        @warn("gammap") = @warn("gamma")
    end
    @warn("gamrat") = if @warn("nst") > 0
            @warn("gamma") / @warn("gammap")
        else
            ONE
        end
end
function CVSetAdams(cv_mem)
    m = @warn("L_MAX")
    M = 3
    local hsum
    if @warn("q") == 1
        (@warn("l"))[1] = ((@warn("l"))[2] = ((@warn("tq"))[2] = ((@warn("tq"))[6] = ONE)))
        (@warn("tq"))[3] = TWO
        (@warn("tq"))[4] = TWELVE
        (@warn("tq"))[5] = CORTES * (@warn("tq"))[3]
        return
    end
    hsum = CVAdamsStart(cv_mem, m)
    M[1] = CVAltSum(@warn("q") - 1, m, 1)
    M[2] = CVAltSum(@warn("q") - 1, m, 2)
    CVAdamsFinish(cv_mem, m, M, hsum)
end
function CVAdamsStart(cv_mem, m)
    local hsum
    local xi_inv
    local sum
    local i
    local j
    hsum = @warn("h")
    m[1] = ONE
    m[i + 1] = ZERO
    begin
        if j == @warn("q") - 1 && @warn("qwait") == 1
            sum = CVAltSum(@warn("q") - 2, m, 2)
            (@warn("tq"))[2] = m[(@warn("q") - 2) + 1] / (@warn("q") * sum)
        end
        xi_inv = @warn("h") / hsum
        m[i + 1] += m[(i - 1) + 1] * xi_inv
        hsum += (@warn("tau"))[j + 1]
    end
    return hsum
end
function CVAdamsFinish(cv_mem, m, M, hsum)
    local i
    local M0_inv
    local xi
    local xi_inv
    M0_inv = ONE / M[1]
    (@warn("l"))[1] = ONE
    (@warn("l"))[i + 1] = M0_inv * (m[(i - 1) + 1] / i)
    xi = hsum / @warn("h")
    xi_inv = ONE / xi
    (@warn("tq"))[3] = (xi * M[1]) / M[2]
    (@warn("tq"))[6] = xi / (@warn("l"))[@warn("q") + 1]
    if @warn("qwait") == 1
        m[i + 1] += m[(i - 1) + 1] * xi_inv
        M[3] = CVAltSum(@warn("q"), m, 2)
        (@warn("tq"))[4] = (@warn("L") * M[1]) / M[3]
    end
    (@warn("tq"))[5] = CORTES * (@warn("tq"))[3]
end
function CVAltSum(iend, a, k)
    local i
    local sign
    local sum
    if iend < 0
        return ZERO
    end
    sum = ZERO
    sign = 1
    begin
        sum += sign * (a[i + 1] / (i + k))
        sign = -sign
    end
    return sum
end
function CVSetBDF(cv_mem)
    local alpha0
    local alpha0_hat
    local xi_inv
    local xistar_inv
    local hsum
    local i
    local j
    (@warn("l"))[1] = ((@warn("l"))[2] = (xi_inv = (xistar_inv = ONE)))
    (@warn("l"))[i + 1] = ZERO
    alpha0 = (alpha0_hat = -ONE)
    hsum = @warn("h")
    if @warn("q") > 1
        begin
            hsum += (@warn("tau"))[(j - 1) + 1]
            xi_inv = @warn("h") / hsum
            alpha0 -= ONE / j
            (@warn("l"))[i + 1] += (@warn("l"))[(i - 1) + 1] * xi_inv
        end
        alpha0 -= ONE / @warn("q")
        xistar_inv = -((@warn("l"))[2]) - alpha0
        hsum += (@warn("tau"))[(@warn("q") - 1) + 1]
        xi_inv = @warn("h") / hsum
        alpha0_hat = -((@warn("l"))[2]) - xi_inv
        (@warn("l"))[i + 1] += (@warn("l"))[(i - 1) + 1] * xistar_inv
    end
    CVSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv)
end
function CVSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv)
    local A1
    local A2
    local A3
    local A4
    local A5
    local A6
    local C
    local CPrime
    local CPrimePrime
    A1 = (ONE - alpha0_hat) + alpha0
    A2 = ONE + @warn("q") * A1
    (@warn("tq"))[3] = @warn("ABS(alpha0*(A2/A1))")
    (@warn("tq"))[6] = @warn("ABS((A2)/(l[q]*xi_inv/xistar_inv))")
    if @warn("qwait") == 1
        C = xistar_inv / (@warn("l"))[@warn("q") + 1]
        A3 = alpha0 + ONE / @warn("q")
        A4 = alpha0_hat + xi_inv
        CPrime = A3 / ((ONE - A4) + A3)
        (@warn("tq"))[2] = @warn("ABS(CPrime/C)")
        hsum += (@warn("tau"))[@warn("q") + 1]
        xi_inv = @warn("h") / hsum
        A5 = alpha0 - ONE / (@warn("q") + 1)
        A6 = alpha0_hat - xi_inv
        CPrimePrime = A2 / ((ONE - A6) + A5)
        (@warn("tq"))[4] = @warn("ABS(CPrimePrime*xi_inv*(q+2)*A5)")
    end
    (@warn("tq"))[5] = CORTES * (@warn("tq"))[3]
end
function CVnls(cv_mem, nflag)
    @switch @warn("iter") begin
            @case FUNCTIONAL
            return CVnlsFunctional(cv_mem)
            @case NEWTON
            return CVnlsNewton(cv_mem, nflag)
        end
end
function CVnlsFunctional(cv_mem)
    local m
    local del
    local delp
    local dcon
    @warn("crate") = ONE
    m = 0
    (@warn("N"), @warn("tn"), (@warn("zn"))[1], @warn("tempv"), @warn("f_data"))
    @warn("nfe") += 1
    N_VConst(ZERO, @warn("acor"))
    begin
        N_VLinearSum(@warn("h"), @warn("tempv"), -ONE, (@warn("zn"))[2], @warn("tempv"))
        N_VScale(@warn("rl1"), @warn("tempv"), @warn("tempv"))
        N_VLinearSum(ONE, (@warn("zn"))[1], ONE, @warn("tempv"), @warn("y"))
        N_VLinearSum(ONE, @warn("tempv"), -ONE, @warn("acor"), @warn("acor"))
        del = N_VWrmsNorm(@warn("acor"), @warn("ewt"))
        N_VScale(ONE, @warn("tempv"), @warn("acor"))
        if m > 0
            @warn("crate") = @warn("MAX(CRDOWN*crate,del/delp)")
        end
        dcon = (del * @warn("MIN(ONE,crate)")) / (@warn("tq"))[5]
        if dcon <= ONE
            @warn("acnrm") = if m == 0
                    del
                else
                    N_VWrmsNorm(@warn("acor"), @warn("ewt"))
                end
            return SOLVED
        end
        m += 1
        if m == @warn("maxcor") || m >= 2 && del > RDIV * delp
            return ?_?(CONV_FAIL)
        end
        delp = del
        (@warn("N"), @warn("tn"), @warn("y"), @warn("tempv"), @warn("f_data"))
        @warn("nfe") += 1
    end
end
function CVnlsNewton(cv_mem, nflag)
    local vtemp1
    local vtemp2
    local vtemp3
    local convfail
    local ier
    local callSetup
    vtemp1 = @warn("acor")
    vtemp2 = @warn("y")
    vtemp3 = @warn("tempv")
    convfail = if nflag == FIRST_CALL || nflag == ?_?(PREV_ERR_FAIL)
            NO_FAILURES
        else
            FAIL_OTHER
        end
    if @warn("setupNonNull")
        callSetup = (((nflag == ?_?(PREV_CONV_FAIL) || nflag == ?_?(PREV_ERR_FAIL)) || @warn("nst") == 0) || @warn("nst") >= @warn("nstlp") + MSBP) || @warn("ABS(gamrat-ONE)") > DGMAX
    else
        @warn("crate") = ONE
        callSetup = FALSE
    end
    begin
        (@warn("N"), @warn("tn"), (@warn("zn"))[1], @warn("ftemp"), @warn("f_data"))
        @warn("nfe") += 1
        if callSetup
            ier = @c((cv_mem, convfail, (@warn("zn"))[1], @warn("ftemp"), &(@warn("jcur")), vtemp1, vtemp2, vtemp3))
            @warn("nsetups") += 1
            callSetup = FALSE
            @warn("gamrat") = (@warn("crate") = ONE)
            @warn("gammap") = @warn("gamma")
            @warn("nstlp") = @warn("nst")
            if ier < 0
                return ?_?(SETUP_FAIL_UNREC)
            end
            if ier > 0
                return ?_?(CONV_FAIL)
            end
        end
        N_VConst(ZERO, @warn("acor"))
        N_VScale(ONE, (@warn("zn"))[1], @warn("y"))
        ier = CVNewtonIteration(cv_mem)
        if ier != TRY_AGAIN
            return ier
        end
        callSetup = TRUE
        convfail = FAIL_BAD_J
    end
end
function CVNewtonIteration(cv_mem)
    local m
    local ret
    local del
    local delp
    local dcon
    local b
    @warn("mnewt") = (m = 0)
    begin
        N_VLinearSum(@warn("rl1"), (@warn("zn"))[2], ONE, @warn("acor"), @warn("tempv"))
        N_VLinearSum(@warn("gamma"), @warn("ftemp"), -ONE, @warn("tempv"), @warn("tempv"))
        b = @warn("tempv")
        ret = (cv_mem, b, @warn("y"), @warn("ftemp"))
        @warn("nni") += 1
        if ret < 0
            return ?_?(SOLVE_FAIL_UNREC)
        end
        if ret > 0
            if !(@warn("jcur")) && @warn("setupNonNull")
                return TRY_AGAIN
            end
            return ?_?(CONV_FAIL)
        end
        del = N_VWrmsNorm(b, @warn("ewt"))
        N_VLinearSum(ONE, @warn("acor"), ONE, b, @warn("acor"))
        N_VLinearSum(ONE, (@warn("zn"))[1], ONE, @warn("acor"), @warn("y"))
        if m > 0
            @warn("crate") = @warn("MAX(CRDOWN*crate,del/delp)")
        end
        dcon = (del * @warn("MIN(ONE,crate)")) / (@warn("tq"))[5]
        if dcon <= ONE
            @warn("acnrm") = if m == 0
                    del
                else
                    N_VWrmsNorm(@warn("acor"), @warn("ewt"))
                end
            @warn("jcur") = FALSE
            return SOLVED
        end
        @warn("mnewt") = @++(m)
        if m == @warn("maxcor") || m >= 2 && del > RDIV * delp
            if !(@warn("jcur")) && @warn("setupNonNull")
                return TRY_AGAIN
            end
            return ?_?(CONV_FAIL)
        end
        delp = del
        (@warn("N"), @warn("tn"), @warn("y"), @warn("ftemp"), @warn("f_data"))
        @warn("nfe") += 1
    end
end
function CVHandleNFlag(cv_mem, nflagPtr, saved_t, ncfPtr)
    local nflag
    nflag = nflagPtr[]
    if nflag == SOLVED
        return DO_ERROR_TEST
    end
    @warn("ncfn") += 1
    CVRestore(cv_mem, saved_t)
    if nflag == ?_?(SETUP_FAIL_UNREC)
        return ?_?(SETUP_FAILED)
    end
    if nflag == ?_?(SOLVE_FAIL_UNREC)
        return ?_?(SOLVE_FAILED)
    end
    (ncfPtr[])[]
    @warn("etamax") = ONE
    if @warn("ABS(h)") <= @warn("hmin") * ONEPSM || ncfPtr[] == MXNCF
        return ?_?(REP_CONV_FAIL)
    end
    @warn("eta") = @warn("MAX(ETACF,hmin/ABS(h))")
    nflagPtr[] = ?_?(PREV_CONV_FAIL)
    CVRescale(cv_mem)
    return ?_?(PREDICT_AGAIN)
end
function CVRestore(cv_mem, saved_t)
    local j
    local k
    @warn("tn") = saved_t
    N_VLinearSum(ONE, (@warn("zn"))[(j - 1) + 1], -ONE, (@warn("zn"))[j + 1], (@warn("zn"))[(j - 1) + 1])
end
function CVDoErrorTest(cv_mem, nflagPtr, kflagPtr, saved_t, nefPtr, dsmPtr)
    local dsm
    dsm = @warn("acnrm") / (@warn("tq"))[3]
    dsmPtr[] = dsm
    if dsm <= ONE
        return TRUE
    end
    (nefPtr[])[]
    @warn("netf") += 1
    nflagPtr[] = ?_?(PREV_ERR_FAIL)
    CVRestore(cv_mem, saved_t)
    if @warn("ABS(h)") <= @warn("hmin") * ONEPSM || nefPtr[] == MXNEF
        kflagPtr[] = ?_?(REP_ERR_FAIL)
        return FALSE
    end
    @warn("etamax") = ONE
    if nefPtr[] <= MXNEF1
        @warn("eta") = ONE / (RPowerR(BIAS2 * dsm, ONE / @warn("L")) + ADDON)
        @warn("eta") = @warn("MAX(ETAMIN,MAX(eta,hmin/ABS(h)))")
        if nefPtr[] >= SMALL_NEF
            @warn("eta") = @warn("MIN(eta,ETAMXF)")
        end
        CVRescale(cv_mem)
        return FALSE
    end
    if @warn("q") > 1
        @warn("eta") = @warn("MAX(ETAMIN,hmin/ABS(h))")
        CVAdjustOrder(cv_mem, -1)
        @warn("L") = @warn("q")
        @warn("q") -= 1
        @warn("qwait") = @warn("L")
        CVRescale(cv_mem)
        return FALSE
    end
    @warn("eta") = @warn("MAX(ETAMIN,hmin/ABS(h))")
    @warn("h") *= @warn("eta")
    @warn("hscale") = @warn("h")
    @warn("qwait") = LONG_WAIT
    (@warn("N"), @warn("tn"), (@warn("zn"))[1], @warn("tempv"), @warn("f_data"))
    @warn("nfe") += 1
    N_VScale(@warn("h"), @warn("tempv"), (@warn("zn"))[2])
    return FALSE
end
function CVCompleteStep(cv_mem)
    local i
    local j
    @warn("nst") += 1
    @warn("hu") = @warn("h")
    @warn("qu") = @warn("q")
    (@warn("tau"))[i + 1] = (@warn("tau"))[(i - 1) + 1]
    if @warn("q") == 1 && @warn("nst") > 1
        (@warn("tau"))[3] = (@warn("tau"))[2]
    end
    (@warn("tau"))[2] = @warn("h")
    N_VLinearSum((@warn("l"))[j + 1], @warn("acor"), ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
    @warn("qwait") -= 1
    if @warn("qwait") == 1 && @warn("q") != @warn("qmax")
        N_VScale(ONE, @warn("acor"), (@warn("zn"))[@warn("qmax") + 1])
        @warn("saved_tq5") = (@warn("tq"))[6]
    end
end
function CVPrepareNextStep(cv_mem, dsm)
    local etaqm1
    local etaq
    local etaqp1
    if @warn("etamax") == ONE
        @warn("qwait") = @warn("MAX(qwait,2)")
        @warn("qprime") = @warn("q")
        @warn("hprime") = @warn("h")
        @warn("eta") = ONE
        @warn("etamax") = if @warn("nst") <= SMALL_NST
                ETAMX2
            else
                ETAMX3
            end
        N_VScale(ONE / (@warn("tq"))[3], @warn("acor"), @warn("acor"))
        return
    end
    etaq = ONE / (RPowerR(BIAS2 * dsm, ONE / @warn("L")) + ADDON)
    if @warn("qwait") != 0
        @warn("eta") = etaq
        @warn("qprime") = @warn("q")
        CVSetEta(cv_mem)
        return
    end
    @warn("qwait") = 2
    etaqm1 = CVComputeEtaqm1(cv_mem)
    etaqp1 = CVComputeEtaqp1(cv_mem)
    CVChooseEta(cv_mem, etaqm1, etaq, etaqp1)
    CVSetEta(cv_mem)
end
function CVSetEta(cv_mem)
    if @warn("eta") < THRESH
        @warn("eta") = ONE
        @warn("hprime") = @warn("h")
    else
        @warn("eta") = @warn("MIN(eta,etamax)")
        @warn("eta") /= @warn("MAX(ONE,ABS(h)*hmax_inv*eta)")
        @warn("hprime") = @warn("h") * @warn("eta")
    end
    @warn("etamax") = if @warn("nst") <= SMALL_NST
            ETAMX2
        else
            ETAMX3
        end
    N_VScale(ONE / (@warn("tq"))[3], @warn("acor"), @warn("acor"))
end
function CVComputeEtaqm1(cv_mem)
    local etaqm1
    local ddn
    etaqm1 = ZERO
    if @warn("q") > 1
        ddn = N_VWrmsNorm((@warn("zn"))[@warn("q") + 1], @warn("ewt")) / (@warn("tq"))[2]
        etaqm1 = ONE / (RPowerR(BIAS1 * ddn, ONE / @warn("q")) + ADDON)
    end
    return etaqm1
end
function CVComputeEtaqp1(cv_mem)
    local etaqp1
    local dup
    local cquot
    etaqp1 = ZERO
    if @warn("q") != @warn("qmax")
        cquot = ((@warn("tq"))[6] / @warn("saved_tq5")) * RPowerI(@warn("h") / (@warn("tau"))[3], @warn("L"))
        N_VLinearSum(-cquot, (@warn("zn"))[@warn("qmax") + 1], ONE, @warn("acor"), @warn("tempv"))
        dup = N_VWrmsNorm(@warn("tempv"), @warn("ewt")) / (@warn("tq"))[4]
        etaqp1 = ONE / (RPowerR(BIAS3 * dup, ONE / (@warn("L") + 1)) + ADDON)
    end
    return etaqp1
end
function CVChooseEta(cv_mem, etaqm1, etaq, etaqp1)
    local etam
    etam = @warn("MAX(etaqm1,MAX(etaq,etaqp1))")
    if etam < THRESH
        @warn("eta") = ONE
        @warn("qprime") = @warn("q")
        return
    end
    if etam == etaq
        @warn("eta") = etaq
        @warn("qprime") = @warn("q")
    elseif etam == etaqm1
        @warn("eta") = etaqm1
        @warn("qprime") = @warn("q") - 1
    else
        @warn("eta") = etaqp1
        @warn("qprime") = @warn("q") + 1
        N_VScale(ONE, @warn("acor"), (@warn("zn"))[@warn("qmax") + 1])
    end
end
function CVHandleFailure(cv_mem, kflag)
    N_VProd(@warn("acor"), @warn("ewt"), @warn("tempv"))
    N_VAbs(@warn("tempv"), @warn("tempv"))
    @switch kflag begin
            @case ?_?(REP_ERR_FAIL)
            fprintf(@warn("errfp"), MSG_ERR_FAILS, @warn("tn"), @warn("h"))
            return ERR_FAILURE
            @case ?_?(REP_CONV_FAIL)
            fprintf(@warn("errfp"), MSG_CONV_FAILS, @warn("tn"), @warn("h"))
            return CONV_FAILURE
            @case ?_?(SETUP_FAILED)
            fprintf(@warn("errfp"), MSG_SETUP_FAILED, @warn("tn"))
            return SETUP_FAILURE
            @case ?_?(SOLVE_FAILED)
            fprintf(@warn("errfp"), MSG_SOLVE_FAILED, @warn("tn"))
            return SOLVE_FAILURE
        end
end
function CVodeMalloc(N, f, t0, y0, lmm, iter, itol, reltol, abstol, f_data, errfp, optIn, iopt, ropt, machEnv)
    local allocOK
    local ioptExists
    local roptExists
    local neg_abstol
    local ewtsetOK
    local maxord
    local cv_mem
    local fp
    fp = if errfp == @warn("NULL")
            __stdoutp
        else
            errfp
        end
    if y0 == @warn("NULL")
        fprintf(fp, MSG_Y0_NULL)
        return @warn("NULL")
    end
    if N <= 0
        fprintf(fp, MSG_BAD_N, N)
        return @warn("NULL")
    end
    if lmm != ADAMS && lmm != BDF
        fprintf(fp, MSG_BAD_LMM, lmm, ADAMS, BDF)
        return @warn("NULL")
    end
    if iter != FUNCTIONAL && iter != NEWTON
        fprintf(fp, MSG_BAD_ITER, iter, FUNCTIONAL, NEWTON)
        return @warn("NULL")
    end
    if itol != SS && itol != SV
        fprintf(fp, MSG_BAD_ITOL, itol, SS, SV)
        return @warn("NULL")
    end
    if f == @warn("NULL")
        fprintf(fp, MSG_F_NULL)
        return @warn("NULL")
    end
    if reltol == @warn("NULL")
        fprintf(fp, MSG_RELTOL_NULL)
        return @warn("NULL")
    end
    if reltol[] < ZERO
        fprintf(fp, MSG_BAD_RELTOL, reltol[])
        return @warn("NULL")
    end
    if abstol == @warn("NULL")
        fprintf(fp, MSG_ABSTOL_NULL)
        return @warn("NULL")
    end
    if itol == SS
        neg_abstol = (Ptr{real}(abstol))[] < ZERO
    else
        neg_abstol = N_VMin(N_Vector(abstol)) < ZERO
    end
    if neg_abstol
        fprintf(fp, MSG_BAD_ABSTOL)
        return @warn("NULL")
    end
    if optIn != FALSE && optIn != TRUE
        fprintf(fp, MSG_BAD_OPTIN, optIn, FALSE, TRUE)
        return @warn("NULL")
    end
    if (optIn && iopt == @warn("NULL")) && ropt == @warn("NULL")
        fprintf(fp, MSG_BAD_OPT)
        return @warn("NULL")
    end
    ioptExists = iopt != @warn("NULL")
    roptExists = ropt != @warn("NULL")
    if optIn && roptExists
        if ropt[HMAX + 1] > ZERO && ropt[HMIN + 1] > ropt[HMAX + 1]
            fprintf(fp, MSG_BAD_HMIN_HMAX, ropt[HMIN + 1], ropt[HMAX + 1])
            return @warn("NULL")
        end
    end
    maxord = if lmm == ADAMS
            ADAMS_Q_MAX
        else
            BDF_Q_MAX
        end
    if optIn && ioptExists
        if iopt[MAXORD + 1] > 0
            maxord = @warn("MIN(maxord,iopt[MAXORD])")
        end
    end
    cv_mem = CVodeMem(malloc(CVodeMemRec()))
    if cv_mem == @warn("NULL")
        fprintf(fp, MSG_MEM_FAIL)
        return @warn("NULL")
    end
    allocOK = CVAllocVectors(cv_mem, N, maxord, machEnv)
    if !allocOK
        fprintf(fp, MSG_MEM_FAIL)
        free(cv_mem)
        return @warn("NULL")
    end
    ewtsetOK = CVEwtSet(cv_mem, reltol, abstol, itol, y0, N)
    if !ewtsetOK
        fprintf(fp, MSG_BAD_EWT)
        CVFreeVectors(cv_mem, maxord)
        free(cv_mem)
        return @warn("NULL")
    end
    cv_mem.cv_N = N
    cv_mem.cv_f = f
    cv_mem.cv_f_data = f_data
    cv_mem.cv_lmm = lmm
    cv_mem.cv_iter = iter
    cv_mem.cv_itol = itol
    cv_mem.cv_reltol = reltol
    cv_mem.cv_abstol = abstol
    cv_mem.cv_iopt = iopt
    cv_mem.cv_ropt = ropt
    cv_mem.cv_errfp = fp
    @warn("tn") = t0
    @warn("machenv") = machEnv
    @warn("q") = 1
    @warn("L") = 2
    @warn("qwait") = @warn("L")
    @warn("qmax") = maxord
    @warn("etamax") = ETAMX1
    @warn("uround") = UnitRoundoff()
    @warn("linit") = @warn("NULL")
    @warn("lsetup") = @warn("NULL")
    @warn("lsolve") = @warn("NULL")
    @warn("lfree") = @warn("NULL")
    @warn("lmem") = @warn("NULL")
    @warn("linitOK") = FALSE
    N_VScale(ONE, y0, (@warn("zn"))[1])
    @warn("hmin") = HMIN_DEFAULT
    @warn("hmax_inv") = HMAX_INV_DEFAULT
    if optIn && roptExists
        if ropt[HMIN + 1] > ZERO
            @warn("hmin") = ropt[HMIN + 1]
        end
        if ropt[HMAX + 1] > ZERO
            @warn("hmax_inv") = ONE / ropt[HMAX + 1]
        end
    end
    @warn("mxhnil") = MXHNIL_DEFAULT
    @warn("mxstep") = MXSTEP_DEFAULT
    if optIn && ioptExists
        if iopt[MXHNIL + 1] > 0
            @warn("mxhnil") = iopt[MXHNIL + 1]
        end
        if iopt[MXSTEP + 1] > 0
            @warn("mxstep") = iopt[MXSTEP + 1]
        end
    end
    if !optIn && roptExists
        ropt[H0 + 1] = ZERO
    end
    @warn("maxcor") = if iter == NEWTON
            NEWT_MAXCOR
        else
            FUNC_MAXCOR
        end
    @warn("nst") = (@warn("nfe") = (@warn("ncfn") = (@warn("netf") = (@warn("nni") = (@warn("nsetups") = (@warn("nhnil") = (@warn("nstlp") = 0)))))))
    @warn("qu") = 0
    @warn("hu") = ZERO
    @warn("tolsf") = ONE
    if ioptExists
        iopt[NST + 1] = (iopt[NFE + 1] = (iopt[NSETUPS + 1] = (iopt[NNI + 1] = 0)))
        iopt[NCFN + 1] = (iopt[NETF + 1] = 0)
        iopt[QU + 1] = @warn("qu")
        iopt[QCUR + 1] = 0
        iopt[LENRW + 1] = @warn("lrw")
        iopt[LENIW + 1] = @warn("liw")
    end
    if roptExists
        ropt[HU + 1] = @warn("hu")
        ropt[HCUR + 1] = ZERO
        ropt[TCUR + 1] = t0
        ropt[TOLSF + 1] = @warn("tolsf")
    end
    return Ptr{Cvoid}(cv_mem)
end
function CVReInit(cvode_mem, f, t0, y0, lmm, iter, itol, reltol, abstol, f_data, errfp, optIn, iopt, ropt, machEnv)
    local allocOK
    local ioptExists
    local roptExists
    local neg_abstol
    local ewtsetOK
    local maxord
    local cv_mem
    local fp
    local N
    fp = if errfp == @warn("NULL")
            __stdoutp
        else
            errfp
        end
    if cvode_mem == @warn("NULL")
        fprintf(fp, MSG_REI_NO_MEM)
        return CVREI_NO_MEM
    end
    cv_mem = CVodeMem(cvode_mem)
    if y0 == @warn("NULL")
        fprintf(fp, MSG_Y0_NULL)
        return CVREI_ILL_INPUT
    end
    if lmm != ADAMS && lmm != BDF
        fprintf(fp, MSG_BAD_LMM, lmm, ADAMS, BDF)
        return CVREI_ILL_INPUT
    end
    if iter != FUNCTIONAL && iter != NEWTON
        fprintf(fp, MSG_BAD_ITER, iter, FUNCTIONAL, NEWTON)
        return CVREI_ILL_INPUT
    end
    if itol != SS && itol != SV
        fprintf(fp, MSG_BAD_ITOL, itol, SS, SV)
        return CVREI_ILL_INPUT
    end
    if f == @warn("NULL")
        fprintf(fp, MSG_F_NULL)
        return CVREI_ILL_INPUT
    end
    if reltol == @warn("NULL")
        fprintf(fp, MSG_RELTOL_NULL)
        return CVREI_ILL_INPUT
    end
    if reltol[] < ZERO
        fprintf(fp, MSG_BAD_RELTOL, reltol[])
        return CVREI_ILL_INPUT
    end
    if abstol == @warn("NULL")
        fprintf(fp, MSG_ABSTOL_NULL)
        return CVREI_ILL_INPUT
    end
    if itol == SS
        neg_abstol = (Ptr{real}(abstol))[] < ZERO
    else
        neg_abstol = N_VMin(N_Vector(abstol)) < ZERO
    end
    if neg_abstol
        fprintf(fp, MSG_BAD_ABSTOL)
        return CVREI_ILL_INPUT
    end
    if optIn != FALSE && optIn != TRUE
        fprintf(fp, MSG_BAD_OPTIN, optIn, FALSE, TRUE)
        return CVREI_ILL_INPUT
    end
    if (optIn && iopt == @warn("NULL")) && ropt == @warn("NULL")
        fprintf(fp, MSG_BAD_OPT)
        return CVREI_ILL_INPUT
    end
    ioptExists = iopt != @warn("NULL")
    roptExists = ropt != @warn("NULL")
    if optIn && roptExists
        if ropt[HMAX + 1] > ZERO && ropt[HMIN + 1] > ropt[HMAX + 1]
            fprintf(fp, MSG_BAD_HMIN_HMAX, ropt[HMIN + 1], ropt[HMAX + 1])
            return CVREI_ILL_INPUT
        end
    end
    maxord = if lmm == ADAMS
            ADAMS_Q_MAX
        else
            BDF_Q_MAX
        end
    if optIn && ioptExists
        if iopt[MAXORD + 1] > 0
            maxord = @warn("MIN(maxord,iopt[MAXORD])")
        end
    end
    if maxord > @warn("qmax")
        fprintf(fp, MSG_REI_MAXORD, @warn("qmax"), maxord)
        return CVREI_ILL_INPUT
    end
    N = cv_mem.cv_N
    ewtsetOK = CVEwtSet(cv_mem, reltol, abstol, itol, y0, N)
    if !ewtsetOK
        fprintf(fp, MSG_BAD_EWT)
        return CVREI_ILL_INPUT
    end
    cv_mem.cv_f = f
    cv_mem.cv_f_data = f_data
    cv_mem.cv_lmm = lmm
    cv_mem.cv_iter = iter
    cv_mem.cv_itol = itol
    cv_mem.cv_reltol = reltol
    cv_mem.cv_abstol = abstol
    cv_mem.cv_iopt = iopt
    cv_mem.cv_ropt = ropt
    cv_mem.cv_errfp = fp
    @warn("tn") = t0
    @warn("machenv") = machEnv
    @warn("q") = 1
    @warn("L") = 2
    @warn("qwait") = @warn("L")
    @warn("qmax") = maxord
    @warn("etamax") = ETAMX1
    @warn("uround") = UnitRoundoff()
    @warn("linit") = @warn("NULL")
    @warn("lsetup") = @warn("NULL")
    @warn("lsolve") = @warn("NULL")
    @warn("lfree") = @warn("NULL")
    @warn("lmem") = @warn("NULL")
    @warn("linitOK") = FALSE
    N_VScale(ONE, y0, (@warn("zn"))[1])
    @warn("hmin") = HMIN_DEFAULT
    @warn("hmax_inv") = HMAX_INV_DEFAULT
    if optIn && roptExists
        if ropt[HMIN + 1] > ZERO
            @warn("hmin") = ropt[HMIN + 1]
        end
        if ropt[HMAX + 1] > ZERO
            @warn("hmax_inv") = ONE / ropt[HMAX + 1]
        end
    end
    @warn("mxhnil") = MXHNIL_DEFAULT
    @warn("mxstep") = MXSTEP_DEFAULT
    if optIn && ioptExists
        if iopt[MXHNIL + 1] > 0
            @warn("mxhnil") = iopt[MXHNIL + 1]
        end
        if iopt[MXSTEP + 1] > 0
            @warn("mxstep") = iopt[MXSTEP + 1]
        end
    end
    if !optIn && roptExists
        ropt[H0 + 1] = ZERO
    end
    @warn("maxcor") = if iter == NEWTON
            NEWT_MAXCOR
        else
            FUNC_MAXCOR
        end
    @warn("nst") = (@warn("nfe") = (@warn("ncfn") = (@warn("netf") = (@warn("nni") = (@warn("nsetups") = (@warn("nhnil") = (@warn("nstlp") = 0)))))))
    @warn("qu") = 0
    @warn("hu") = ZERO
    @warn("tolsf") = ONE
    if ioptExists
        iopt[NST + 1] = (iopt[NFE + 1] = (iopt[NSETUPS + 1] = (iopt[NNI + 1] = 0)))
        iopt[NCFN + 1] = (iopt[NETF + 1] = 0)
        iopt[QU + 1] = @warn("qu")
        iopt[QCUR + 1] = 0
        iopt[LENRW + 1] = @warn("lrw")
        iopt[LENIW + 1] = @warn("liw")
    end
    if roptExists
        ropt[HU + 1] = @warn("hu")
        ropt[HCUR + 1] = ZERO
        ropt[TCUR + 1] = t0
        ropt[TOLSF + 1] = @warn("tolsf")
    end
    return SUCCESS
end
function CVode(cvode_mem, tout, yout, t, itask)
    local nstloc
    local kflag
    local istate
    local next_q
    local ier
    local rh
    local next_h
    local hOK
    local ewtsetOK
    local cv_mem
    cv_mem = CVodeMem(cvode_mem)
    if cvode_mem == @warn("NULL")
        fprintf(__stdoutp, MSG_CVODE_NO_MEM)
        return CVODE_NO_MEM
    end
    if (@warn("y") = yout) == @warn("NULL")
        fprintf(@warn("errfp"), MSG_YOUT_NULL)
        return ILL_INPUT
    end
    if t == @warn("NULL")
        fprintf(@warn("errfp"), MSG_T_NULL)
        return ILL_INPUT
    end
    t[] = @warn("tn")
    if itask != NORMAL && itask != ONE_STEP
        fprintf(@warn("errfp"), MSG_BAD_ITASK, itask, NORMAL, ONE_STEP)
        return ILL_INPUT
    end
    if @warn("nst") == 0
        if @warn("iter") == NEWTON
            if @warn("linit") == @warn("NULL")
                fprintf(@warn("errfp"), MSG_LINIT_NULL)
                return ILL_INPUT
            end
            if @warn("lsetup") == @warn("NULL")
                fprintf(@warn("errfp"), MSG_LSETUP_NULL)
                return ILL_INPUT
            end
            if @warn("lsolve") == @warn("NULL")
                fprintf(@warn("errfp"), MSG_LSOLVE_NULL)
                return ILL_INPUT
            end
            if @warn("lfree") == @warn("NULL")
                fprintf(@warn("errfp"), MSG_LFREE_NULL)
                return ILL_INPUT
            end
            @warn("linitOK") = @c((cv_mem, &(@warn("setupNonNull")))) == LINIT_OK
            if !(@warn("linitOK"))
                fprintf(@warn("errfp"), MSG_LINIT_FAIL)
                return ILL_INPUT
            end
        end
        (@warn("N"), @warn("tn"), (@warn("zn"))[1], (@warn("zn"))[2], @warn("f_data"))
        @warn("nfe") = 1
        @warn("h") = ZERO
        if @warn("ropt") != @warn("NULL")
            @warn("h") = (@warn("ropt"))[H0 + 1]
        end
        if @warn("h") != ZERO && (tout - @warn("tn")) * @warn("h") < ZERO
            fprintf(@warn("errfp"), MSG_BAD_H0, @warn("h"), tout - @warn("tn"))
            return ILL_INPUT
        end
        if @warn("h") == ZERO
            hOK = CVHin(cv_mem, tout)
            if !hOK
                fprintf(@warn("errfp"), MSG_TOO_CLOSE, tout, @warn("tn"))
                return ILL_INPUT
            end
        end
        rh = @warn("ABS(h)") * @warn("hmax_inv")
        if rh > ONE
            @warn("h") /= rh
        end
        if @warn("ABS(h)") < @warn("hmin")
            @warn("h") *= @warn("hmin") / @warn("ABS(h)")
        end
        @warn("hscale") = @warn("h")
        N_VScale(@warn("h"), (@warn("zn"))[2], (@warn("zn"))[2])
    end
    if (itask == NORMAL && @warn("nst") > 0) && (@warn("tn") - tout) * @warn("h") >= ZERO
        t[] = tout
        ier = CVodeDky(cv_mem, tout, 0, yout)
        if ier != OKAY
            fprintf(@warn("errfp"), MSG_BAD_TOUT, tout)
            return ILL_INPUT
        end
        return SUCCESS
    end
    nstloc = 0
    begin
        next_h = @warn("h")
        next_q = @warn("q")
        if @warn("nst") > 0
            ewtsetOK = CVEwtSet(cv_mem, @warn("reltol"), @warn("abstol"), @warn("itol"), (@warn("zn"))[1], @warn("N"))
            if !ewtsetOK
                fprintf(@warn("errfp"), MSG_EWT_NOW_BAD, @warn("tn"))
                istate = ILL_INPUT
                t[] = @warn("tn")
                N_VScale(ONE, (@warn("zn"))[1], yout)
                break
            end
        end
        if nstloc >= @warn("mxstep")
            fprintf(@warn("errfp"), MSG_MAX_STEPS, @warn("tn"), @warn("mxstep"), tout)
            istate = TOO_MUCH_WORK
            t[] = @warn("tn")
            N_VScale(ONE, (@warn("zn"))[1], yout)
            break
        end
        if (@warn("tolsf") = @warn("uround") * N_VWrmsNorm((@warn("zn"))[1], @warn("ewt"))) > ONE
            fprintf(@warn("errfp"), MSG_TOO_MUCH_ACC, @warn("tn"))
            istate = TOO_MUCH_ACC
            t[] = @warn("tn")
            N_VScale(ONE, (@warn("zn"))[1], yout)
            @warn("tolsf") *= TWO
            break
        end
        if @warn("tn") + @warn("h") == @warn("tn")
            @warn("nhnil") += 1
            if @warn("nhnil") <= @warn("mxhnil")
                fprintf(@warn("errfp"), MSG_HNIL, @warn("tn"), @warn("h"))
            end
            if @warn("nhnil") == @warn("mxhnil")
                fprintf(@warn("errfp"), MSG_HNIL_DONE, @warn("mxhnil"))
            end
        end
        kflag = CVStep(cv_mem)
        if kflag != SUCCESS_STEP
            istate = CVHandleFailure(cv_mem, kflag)
            t[] = @warn("tn")
            N_VScale(ONE, (@warn("zn"))[1], yout)
            break
        end
        nstloc += 1
        if itask == ONE_STEP
            istate = SUCCESS
            t[] = @warn("tn")
            N_VScale(ONE, (@warn("zn"))[1], yout)
            next_q = @warn("qprime")
            next_h = @warn("hprime")
            break
        end
        if (@warn("tn") - tout) * @warn("h") >= ZERO
            istate = SUCCESS
            t[] = tout
            CVodeDky(cv_mem, tout, 0, yout)
            next_q = @warn("qprime")
            next_h = @warn("hprime")
            break
        end
    end
    if @warn("iopt") != @warn("NULL")
        (@warn("iopt"))[NST + 1] = @warn("nst")
        (@warn("iopt"))[NFE + 1] = @warn("nfe")
        (@warn("iopt"))[NSETUPS + 1] = @warn("nsetups")
        (@warn("iopt"))[NNI + 1] = @warn("nni")
        (@warn("iopt"))[NCFN + 1] = @warn("ncfn")
        (@warn("iopt"))[NETF + 1] = @warn("netf")
        (@warn("iopt"))[QU + 1] = @warn("q")
        (@warn("iopt"))[QCUR + 1] = next_q
    end
    if @warn("ropt") != @warn("NULL")
        (@warn("ropt"))[HU + 1] = @warn("h")
        (@warn("ropt"))[HCUR + 1] = next_h
        (@warn("ropt"))[TCUR + 1] = @warn("tn")
        (@warn("ropt"))[TOLSF + 1] = @warn("tolsf")
    end
    return istate
end
function CVodeDky(cvode_mem, t, k, dky)
    local s
    local c
    local r
    local tfuzz
    local tp
    local tn1
    local i
    local j
    local cv_mem
    cv_mem = CVodeMem(cvode_mem)
    if cvode_mem == @warn("NULL")
        fprintf(__stdoutp, MSG_DKY_NO_MEM)
        return DKY_NO_MEM
    end
    if dky == @warn("NULL")
        fprintf(__stdoutp, MSG_BAD_DKY)
        return BAD_DKY
    end
    if k < 0 || k > @warn("q")
        fprintf(@warn("errfp"), MSG_BAD_K, k)
        return BAD_K
    end
    tfuzz = (FUZZ_FACTOR * @warn("uround")) * (@warn("tn") + @warn("hu"))
    tp = (@warn("tn") - @warn("hu")) - tfuzz
    tn1 = @warn("tn") + tfuzz
    if (t - tp) * (t - tn1) > ZERO
        fprintf(@warn("errfp"), MSG_BAD_T, t, @warn("tn") - @warn("hu"), @warn("tn"))
        return BAD_T
    end
    s = (t - @warn("tn")) / @warn("h")
    begin
        c = ONE
        c *= i
        if j == @warn("q")
            N_VScale(c, (@warn("zn"))[@warn("q") + 1], dky)
        else
            N_VLinearSum(c, (@warn("zn"))[j + 1], s, dky, dky)
        end
    end
    if k == 0
        return OKAY
    end
    r = RPowerI(@warn("h"), -k)
    N_VScale(r, dky, dky)
    return OKAY
end
function CVodeFree(cvode_mem)
    local cv_mem
    cv_mem = CVodeMem(cvode_mem)
    if cvode_mem == @warn("NULL")
        return
    end
    CVFreeVectors(cv_mem, @warn("qmax"))
    if @warn("iter") == NEWTON && @warn("linitOK")
        (cv_mem)
    end
    free(cv_mem)
end
function CVAllocVectors(cv_mem, neq, maxord, machEnv)
    local i
    local j
    @warn("ewt") = N_VNew(neq, machEnv)
    if @warn("ewt") == @warn("NULL")
        return FALSE
    end
    @warn("acor") = N_VNew(neq, machEnv)
    if @warn("acor") == @warn("NULL")
        N_VFree(@warn("ewt"))
        return FALSE
    end
    @warn("tempv") = N_VNew(neq, machEnv)
    if @warn("tempv") == @warn("NULL")
        N_VFree(@warn("ewt"))
        N_VFree(@warn("acor"))
        return FALSE
    end
    @warn("ftemp") = N_VNew(neq, machEnv)
    if @warn("ftemp") == @warn("NULL")
        N_VFree(@warn("tempv"))
        N_VFree(@warn("ewt"))
        N_VFree(@warn("acor"))
        return FALSE
    end
    begin
        (@warn("zn"))[j + 1] = N_VNew(neq, machEnv)
        if (@warn("zn"))[j + 1] == @warn("NULL")
            N_VFree(@warn("ewt"))
            N_VFree(@warn("acor"))
            N_VFree(@warn("tempv"))
            N_VFree(@warn("ftemp"))
            N_VFree((@warn("zn"))[i + 1])
            return FALSE
        end
    end
    @warn("lrw") = (maxord + 5) * neq
    @warn("liw") = 0
    return TRUE
end
function CVFreeVectors(cv_mem, maxord)
    local j
    N_VFree(@warn("ewt"))
    N_VFree(@warn("acor"))
    N_VFree(@warn("tempv"))
    N_VFree(@warn("ftemp"))
    N_VFree((@warn("zn"))[j + 1])
end
function CVEwtSet(cv_mem, rtol, atol, tol_type, ycur, neq)
    @switch tol_type begin
            @case SS
            return CVEwtSetSS(cv_mem, rtol, Ptr{real}(atol), ycur, neq)
            @case SV
            return CVEwtSetSV(cv_mem, rtol, N_Vector(atol), ycur, neq)
        end
end
function CVEwtSetSS(cv_mem, rtol, atol, ycur, neq)
    local rtoli
    local atoli
    rtoli = rtol[]
    atoli = atol[]
    N_VAbs(ycur, @warn("tempv"))
    N_VScale(rtoli, @warn("tempv"), @warn("tempv"))
    N_VAddConst(@warn("tempv"), atoli, @warn("tempv"))
    if N_VMin(@warn("tempv")) <= ZERO
        return FALSE
    end
    N_VInv(@warn("tempv"), @warn("ewt"))
    return TRUE
end
function CVEwtSetSV(cv_mem, rtol, atol, ycur, neq)
    local rtoli
    rtoli = rtol[]
    N_VAbs(ycur, @warn("tempv"))
    N_VLinearSum(rtoli, @warn("tempv"), ONE, atol, @warn("tempv"))
    if N_VMin(@warn("tempv")) <= ZERO
        return FALSE
    end
    N_VInv(@warn("tempv"), @warn("ewt"))
    return TRUE
end
function CVHin(cv_mem, tout)
    local sign
    local count
    local tdiff
    local tdist
    local tround
    local hlb
    local hub
    local hg
    local hgs
    local hnew
    local hrat
    local h0
    local yddnrm
    if (tdiff = tout - @warn("tn")) == ZERO
        return FALSE
    end
    sign = if tdiff > ZERO
            1
        else
            -1
        end
    tdist = @warn("ABS(tdiff)")
    tround = @warn("uround") * @warn("MAX(ABS(tn),ABS(tout))")
    if tdist < TWO * tround
        return FALSE
    end
    hlb = HLB_FACTOR * tround
    hub = CVUpperBoundH0(cv_mem, tdist)
    hg = RSqrt(hlb * hub)
    if hub < hlb
        if sign == -1
            hg = -hg
        end
        @warn("h") = hg
        return TRUE
    end
    count = 0
    begin
        hgs = hg * sign
        yddnrm = CVYddNorm(cv_mem, hgs)
        hnew = if (yddnrm * hub) * hub > TWO
                RSqrt(TWO / yddnrm)
            else
                RSqrt(hg * hub)
            end
        count += 1
        if count >= MAX_ITERS
            break
        end
        hrat = hnew / hg
        if hrat > HALF && hrat < TWO
            break
        end
        if count >= 2 && hrat > TWO
            hnew = hg
            break
        end
        hg = hnew
    end
    h0 = H_BIAS * hnew
    if h0 < hlb
        h0 = hlb
    end
    if h0 > hub
        h0 = hub
    end
    if sign == -1
        h0 = -h0
    end
    @warn("h") = h0
    return TRUE
end
function CVUpperBoundH0(cv_mem, tdist)
    local atoli
    local hub_inv
    local hub
    local vectorAtol
    local temp1
    local temp2
    vectorAtol = @warn("itol") == SV
    if !vectorAtol
        atoli = (Ptr{real}(@warn("abstol")))[]
    end
    temp1 = @warn("tempv")
    temp2 = @warn("acor")
    N_VAbs((@warn("zn"))[1], temp1)
    N_VAbs((@warn("zn"))[2], temp2)
    if vectorAtol
        N_VLinearSum(HUB_FACTOR, temp1, ONE, N_Vector(@warn("abstol")), temp1)
    else
        N_VScale(HUB_FACTOR, temp1, temp1)
        N_VAddConst(temp1, atoli, temp1)
    end
    N_VDiv(temp2, temp1, temp1)
    hub_inv = N_VMaxNorm(temp1)
    hub = HUB_FACTOR * tdist
    if hub * hub_inv > ONE
        hub = ONE / hub_inv
    end
    return hub
end
function CVYddNorm(cv_mem, hg)
    local yddnrm
    N_VLinearSum(hg, (@warn("zn"))[2], ONE, (@warn("zn"))[1], @warn("y"))
    (@warn("N"), @warn("tn") + hg, @warn("y"), @warn("tempv"), @warn("f_data"))
    @warn("nfe") += 1
    N_VLinearSum(ONE, @warn("tempv"), -ONE, (@warn("zn"))[2], @warn("tempv"))
    N_VScale(ONE / hg, @warn("tempv"), @warn("tempv"))
    yddnrm = N_VWrmsNorm(@warn("tempv"), @warn("ewt"))
    return yddnrm
end
function CVStep(cv_mem)
    local saved_t
    local dsm
    local ncf
    local nef
    local nflag
    local kflag
    local passed
    saved_t = @warn("tn")
    ncf = (nef = 0)
    nflag = FIRST_CALL
    if @warn("nst") > 0 && @warn("hprime") != @warn("h")
        CVAdjustParams(cv_mem)
    end
    begin
        CVPredict(cv_mem)
        CVSet(cv_mem)
        nflag = CVnls(cv_mem, nflag)
        kflag = @c(CVHandleNFlag(cv_mem, &nflag, saved_t, &ncf))
        if kflag == ?_?(PREDICT_AGAIN)
            continue
        end
        if kflag != DO_ERROR_TEST
            return kflag
        end
        passed = @c(CVDoErrorTest(cv_mem, &nflag, &kflag, saved_t, &nef, &dsm))
        if !passed && kflag == ?_?(REP_ERR_FAIL)
            return kflag
        end
        if passed
            break
        end
    end
    CVCompleteStep(cv_mem)
    CVPrepareNextStep(cv_mem, dsm)
    return SUCCESS_STEP
end
function CVAdjustParams(cv_mem)
    if @warn("qprime") != @warn("q")
        CVAdjustOrder(cv_mem, @warn("qprime") - @warn("q"))
        @warn("q") = @warn("qprime")
        @warn("L") = @warn("q") + 1
        @warn("qwait") = @warn("L")
    end
    CVRescale(cv_mem)
end
function CVAdjustOrder(cv_mem, deltaq)
    if @warn("q") == 2 && deltaq != 1
        return
    end
    @switch @warn("lmm") begin
            @case ADAMS
            CVAdjustAdams(cv_mem, deltaq)
            break
            @case BDF
            CVAdjustBDF(cv_mem, deltaq)
            break
        end
end
function CVAdjustAdams(cv_mem, deltaq)
    local i
    local j
    local xi
    local hsum
    if deltaq == 1
        N_VConst(ZERO, (@warn("zn"))[@warn("L") + 1])
        return
    end
    (@warn("l"))[i + 1] = ZERO
    (@warn("l"))[2] = ONE
    hsum = ZERO
    begin
        hsum += (@warn("tau"))[j + 1]
        xi = hsum / @warn("hscale")
        (@warn("l"))[i + 1] = (@warn("l"))[i + 1] * xi + (@warn("l"))[(i - 1) + 1]
    end
    (@warn("l"))[(j + 1) + 1] = @warn("q") * ((@warn("l"))[j + 1] / (j + 1))
    N_VLinearSum(-((@warn("l"))[j + 1]), (@warn("zn"))[@warn("q") + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
end
function CVAdjustBDF(cv_mem, deltaq)
    @switch deltaq begin
            @case 1
            CVIncreaseBDF(cv_mem)
            return
            @case -1
            CVDecreaseBDF(cv_mem)
            return
        end
end
function CVIncreaseBDF(cv_mem)
    local alpha0
    local alpha1
    local prod
    local xi
    local xiold
    local hsum
    local A1
    local i
    local j
    (@warn("l"))[i + 1] = ZERO
    (@warn("l"))[3] = (alpha1 = (prod = (xiold = ONE)))
    alpha0 = -ONE
    hsum = @warn("hscale")
    if @warn("q") > 1
        begin
            hsum += (@warn("tau"))[(j + 1) + 1]
            xi = hsum / @warn("hscale")
            prod *= xi
            alpha0 -= ONE / (j + 1)
            alpha1 += ONE / xi
            (@warn("l"))[i + 1] = (@warn("l"))[i + 1] * xiold + (@warn("l"))[(i - 1) + 1]
            xiold = xi
        end
    end
    A1 = (-alpha0 - alpha1) / prod
    N_VScale(A1, (@warn("zn"))[@warn("qmax") + 1], (@warn("zn"))[@warn("L") + 1])
    begin
        N_VLinearSum((@warn("l"))[j + 1], (@warn("zn"))[@warn("L") + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
    end
end
function CVDecreaseBDF(cv_mem)
    local hsum
    local xi
    local i
    local j
    (@warn("l"))[i + 1] = ZERO
    (@warn("l"))[3] = ONE
    hsum = ZERO
    begin
        hsum += (@warn("tau"))[j + 1]
        xi = hsum / @warn("hscale")
        (@warn("l"))[i + 1] = (@warn("l"))[i + 1] * xi + (@warn("l"))[(i - 1) + 1]
    end
    N_VLinearSum(-((@warn("l"))[j + 1]), (@warn("zn"))[@warn("q") + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
end
function CVRescale(cv_mem)
    local j
    local factor
    factor = @warn("eta")
    begin
        N_VScale(factor, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
        factor *= @warn("eta")
    end
    @warn("h") = @warn("hscale") * @warn("eta")
    @warn("hscale") = @warn("h")
end
function CVPredict(cv_mem)
    local j
    local k
    @warn("tn") += @warn("h")
    N_VLinearSum(ONE, (@warn("zn"))[(j - 1) + 1], ONE, (@warn("zn"))[j + 1], (@warn("zn"))[(j - 1) + 1])
end
function CVSet(cv_mem)
    @switch @warn("lmm") begin
            @case ADAMS
            CVSetAdams(cv_mem)
            break
            @case BDF
            CVSetBDF(cv_mem)
            break
        end
    @warn("rl1") = ONE / (@warn("l"))[2]
    @warn("gamma") = @warn("h") * @warn("rl1")
    if @warn("nst") == 0
        @warn("gammap") = @warn("gamma")
    end
    @warn("gamrat") = if @warn("nst") > 0
            @warn("gamma") / @warn("gammap")
        else
            ONE
        end
end
function CVSetAdams(cv_mem)
    m = @warn("L_MAX")
    M = 3
    local hsum
    if @warn("q") == 1
        (@warn("l"))[1] = ((@warn("l"))[2] = ((@warn("tq"))[2] = ((@warn("tq"))[6] = ONE)))
        (@warn("tq"))[3] = TWO
        (@warn("tq"))[4] = TWELVE
        (@warn("tq"))[5] = CORTES * (@warn("tq"))[3]
        return
    end
    hsum = CVAdamsStart(cv_mem, m)
    M[1] = CVAltSum(@warn("q") - 1, m, 1)
    M[2] = CVAltSum(@warn("q") - 1, m, 2)
    CVAdamsFinish(cv_mem, m, M, hsum)
end
function CVAdamsStart(cv_mem, m)
    local hsum
    local xi_inv
    local sum
    local i
    local j
    hsum = @warn("h")
    m[1] = ONE
    m[i + 1] = ZERO
    begin
        if j == @warn("q") - 1 && @warn("qwait") == 1
            sum = CVAltSum(@warn("q") - 2, m, 2)
            (@warn("tq"))[2] = m[(@warn("q") - 2) + 1] / (@warn("q") * sum)
        end
        xi_inv = @warn("h") / hsum
        m[i + 1] += m[(i - 1) + 1] * xi_inv
        hsum += (@warn("tau"))[j + 1]
    end
    return hsum
end
function CVAdamsFinish(cv_mem, m, M, hsum)
    local i
    local M0_inv
    local xi
    local xi_inv
    M0_inv = ONE / M[1]
    (@warn("l"))[1] = ONE
    (@warn("l"))[i + 1] = M0_inv * (m[(i - 1) + 1] / i)
    xi = hsum / @warn("h")
    xi_inv = ONE / xi
    (@warn("tq"))[3] = (xi * M[1]) / M[2]
    (@warn("tq"))[6] = xi / (@warn("l"))[@warn("q") + 1]
    if @warn("qwait") == 1
        m[i + 1] += m[(i - 1) + 1] * xi_inv
        M[3] = CVAltSum(@warn("q"), m, 2)
        (@warn("tq"))[4] = (@warn("L") * M[1]) / M[3]
    end
    (@warn("tq"))[5] = CORTES * (@warn("tq"))[3]
end
function CVAltSum(iend, a, k)
    local i
    local sign
    local sum
    if iend < 0
        return ZERO
    end
    sum = ZERO
    sign = 1
    begin
        sum += sign * (a[i + 1] / (i + k))
        sign = -sign
    end
    return sum
end
function CVSetBDF(cv_mem)
    local alpha0
    local alpha0_hat
    local xi_inv
    local xistar_inv
    local hsum
    local i
    local j
    (@warn("l"))[1] = ((@warn("l"))[2] = (xi_inv = (xistar_inv = ONE)))
    (@warn("l"))[i + 1] = ZERO
    alpha0 = (alpha0_hat = -ONE)
    hsum = @warn("h")
    if @warn("q") > 1
        begin
            hsum += (@warn("tau"))[(j - 1) + 1]
            xi_inv = @warn("h") / hsum
            alpha0 -= ONE / j
            (@warn("l"))[i + 1] += (@warn("l"))[(i - 1) + 1] * xi_inv
        end
        alpha0 -= ONE / @warn("q")
        xistar_inv = -((@warn("l"))[2]) - alpha0
        hsum += (@warn("tau"))[(@warn("q") - 1) + 1]
        xi_inv = @warn("h") / hsum
        alpha0_hat = -((@warn("l"))[2]) - xi_inv
        (@warn("l"))[i + 1] += (@warn("l"))[(i - 1) + 1] * xistar_inv
    end
    CVSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv)
end
function CVSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv)
    local A1
    local A2
    local A3
    local A4
    local A5
    local A6
    local C
    local CPrime
    local CPrimePrime
    A1 = (ONE - alpha0_hat) + alpha0
    A2 = ONE + @warn("q") * A1
    (@warn("tq"))[3] = @warn("ABS(alpha0*(A2/A1))")
    (@warn("tq"))[6] = @warn("ABS((A2)/(l[q]*xi_inv/xistar_inv))")
    if @warn("qwait") == 1
        C = xistar_inv / (@warn("l"))[@warn("q") + 1]
        A3 = alpha0 + ONE / @warn("q")
        A4 = alpha0_hat + xi_inv
        CPrime = A3 / ((ONE - A4) + A3)
        (@warn("tq"))[2] = @warn("ABS(CPrime/C)")
        hsum += (@warn("tau"))[@warn("q") + 1]
        xi_inv = @warn("h") / hsum
        A5 = alpha0 - ONE / (@warn("q") + 1)
        A6 = alpha0_hat - xi_inv
        CPrimePrime = A2 / ((ONE - A6) + A5)
        (@warn("tq"))[4] = @warn("ABS(CPrimePrime*xi_inv*(q+2)*A5)")
    end
    (@warn("tq"))[5] = CORTES * (@warn("tq"))[3]
end
function CVnls(cv_mem, nflag)
    @switch @warn("iter") begin
            @case FUNCTIONAL
            return CVnlsFunctional(cv_mem)
            @case NEWTON
            return CVnlsNewton(cv_mem, nflag)
        end
end
function CVnlsFunctional(cv_mem)
    local m
    local del
    local delp
    local dcon
    @warn("crate") = ONE
    m = 0
    (@warn("N"), @warn("tn"), (@warn("zn"))[1], @warn("tempv"), @warn("f_data"))
    @warn("nfe") += 1
    N_VConst(ZERO, @warn("acor"))
    begin
        N_VLinearSum(@warn("h"), @warn("tempv"), -ONE, (@warn("zn"))[2], @warn("tempv"))
        N_VScale(@warn("rl1"), @warn("tempv"), @warn("tempv"))
        N_VLinearSum(ONE, (@warn("zn"))[1], ONE, @warn("tempv"), @warn("y"))
        N_VLinearSum(ONE, @warn("tempv"), -ONE, @warn("acor"), @warn("acor"))
        del = N_VWrmsNorm(@warn("acor"), @warn("ewt"))
        N_VScale(ONE, @warn("tempv"), @warn("acor"))
        if m > 0
            @warn("crate") = @warn("MAX(CRDOWN*crate,del/delp)")
        end
        dcon = (del * @warn("MIN(ONE,crate)")) / (@warn("tq"))[5]
        if dcon <= ONE
            @warn("acnrm") = if m == 0
                    del
                else
                    N_VWrmsNorm(@warn("acor"), @warn("ewt"))
                end
            return SOLVED
        end
        m += 1
        if m == @warn("maxcor") || m >= 2 && del > RDIV * delp
            return ?_?(CONV_FAIL)
        end
        delp = del
        (@warn("N"), @warn("tn"), @warn("y"), @warn("tempv"), @warn("f_data"))
        @warn("nfe") += 1
    end
end
function CVnlsNewton(cv_mem, nflag)
    local vtemp1
    local vtemp2
    local vtemp3
    local convfail
    local ier
    local callSetup
    vtemp1 = @warn("acor")
    vtemp2 = @warn("y")
    vtemp3 = @warn("tempv")
    convfail = if nflag == FIRST_CALL || nflag == ?_?(PREV_ERR_FAIL)
            NO_FAILURES
        else
            FAIL_OTHER
        end
    if @warn("setupNonNull")
        callSetup = (((nflag == ?_?(PREV_CONV_FAIL) || nflag == ?_?(PREV_ERR_FAIL)) || @warn("nst") == 0) || @warn("nst") >= @warn("nstlp") + MSBP) || @warn("ABS(gamrat-ONE)") > DGMAX
    else
        @warn("crate") = ONE
        callSetup = FALSE
    end
    begin
        (@warn("N"), @warn("tn"), (@warn("zn"))[1], @warn("ftemp"), @warn("f_data"))
        @warn("nfe") += 1
        if callSetup
            ier = @c((cv_mem, convfail, (@warn("zn"))[1], @warn("ftemp"), &(@warn("jcur")), vtemp1, vtemp2, vtemp3))
            @warn("nsetups") += 1
            callSetup = FALSE
            @warn("gamrat") = (@warn("crate") = ONE)
            @warn("gammap") = @warn("gamma")
            @warn("nstlp") = @warn("nst")
            if ier < 0
                return ?_?(SETUP_FAIL_UNREC)
            end
            if ier > 0
                return ?_?(CONV_FAIL)
            end
        end
        N_VConst(ZERO, @warn("acor"))
        N_VScale(ONE, (@warn("zn"))[1], @warn("y"))
        ier = CVNewtonIteration(cv_mem)
        if ier != TRY_AGAIN
            return ier
        end
        callSetup = TRUE
        convfail = FAIL_BAD_J
    end
end
function CVNewtonIteration(cv_mem)
    local m
    local ret
    local del
    local delp
    local dcon
    local b
    @warn("mnewt") = (m = 0)
    begin
        N_VLinearSum(@warn("rl1"), (@warn("zn"))[2], ONE, @warn("acor"), @warn("tempv"))
        N_VLinearSum(@warn("gamma"), @warn("ftemp"), -ONE, @warn("tempv"), @warn("tempv"))
        b = @warn("tempv")
        ret = (cv_mem, b, @warn("y"), @warn("ftemp"))
        @warn("nni") += 1
        if ret < 0
            return ?_?(SOLVE_FAIL_UNREC)
        end
        if ret > 0
            if !(@warn("jcur")) && @warn("setupNonNull")
                return TRY_AGAIN
            end
            return ?_?(CONV_FAIL)
        end
        del = N_VWrmsNorm(b, @warn("ewt"))
        N_VLinearSum(ONE, @warn("acor"), ONE, b, @warn("acor"))
        N_VLinearSum(ONE, (@warn("zn"))[1], ONE, @warn("acor"), @warn("y"))
        if m > 0
            @warn("crate") = @warn("MAX(CRDOWN*crate,del/delp)")
        end
        dcon = (del * @warn("MIN(ONE,crate)")) / (@warn("tq"))[5]
        if dcon <= ONE
            @warn("acnrm") = if m == 0
                    del
                else
                    N_VWrmsNorm(@warn("acor"), @warn("ewt"))
                end
            @warn("jcur") = FALSE
            return SOLVED
        end
        @warn("mnewt") = @++(m)
        if m == @warn("maxcor") || m >= 2 && del > RDIV * delp
            if !(@warn("jcur")) && @warn("setupNonNull")
                return TRY_AGAIN
            end
            return ?_?(CONV_FAIL)
        end
        delp = del
        (@warn("N"), @warn("tn"), @warn("y"), @warn("ftemp"), @warn("f_data"))
        @warn("nfe") += 1
    end
end
function CVHandleNFlag(cv_mem, nflagPtr, saved_t, ncfPtr)
    local nflag
    nflag = nflagPtr[]
    if nflag == SOLVED
        return DO_ERROR_TEST
    end
    @warn("ncfn") += 1
    CVRestore(cv_mem, saved_t)
    if nflag == ?_?(SETUP_FAIL_UNREC)
        return ?_?(SETUP_FAILED)
    end
    if nflag == ?_?(SOLVE_FAIL_UNREC)
        return ?_?(SOLVE_FAILED)
    end
    (ncfPtr[])[]
    @warn("etamax") = ONE
    if @warn("ABS(h)") <= @warn("hmin") * ONEPSM || ncfPtr[] == MXNCF
        return ?_?(REP_CONV_FAIL)
    end
    @warn("eta") = @warn("MAX(ETACF,hmin/ABS(h))")
    nflagPtr[] = ?_?(PREV_CONV_FAIL)
    CVRescale(cv_mem)
    return ?_?(PREDICT_AGAIN)
end
function CVRestore(cv_mem, saved_t)
    local j
    local k
    @warn("tn") = saved_t
    N_VLinearSum(ONE, (@warn("zn"))[(j - 1) + 1], -ONE, (@warn("zn"))[j + 1], (@warn("zn"))[(j - 1) + 1])
end
function CVDoErrorTest(cv_mem, nflagPtr, kflagPtr, saved_t, nefPtr, dsmPtr)
    local dsm
    dsm = @warn("acnrm") / (@warn("tq"))[3]
    dsmPtr[] = dsm
    if dsm <= ONE
        return TRUE
    end
    (nefPtr[])[]
    @warn("netf") += 1
    nflagPtr[] = ?_?(PREV_ERR_FAIL)
    CVRestore(cv_mem, saved_t)
    if @warn("ABS(h)") <= @warn("hmin") * ONEPSM || nefPtr[] == MXNEF
        kflagPtr[] = ?_?(REP_ERR_FAIL)
        return FALSE
    end
    @warn("etamax") = ONE
    if nefPtr[] <= MXNEF1
        @warn("eta") = ONE / (RPowerR(BIAS2 * dsm, ONE / @warn("L")) + ADDON)
        @warn("eta") = @warn("MAX(ETAMIN,MAX(eta,hmin/ABS(h)))")
        if nefPtr[] >= SMALL_NEF
            @warn("eta") = @warn("MIN(eta,ETAMXF)")
        end
        CVRescale(cv_mem)
        return FALSE
    end
    if @warn("q") > 1
        @warn("eta") = @warn("MAX(ETAMIN,hmin/ABS(h))")
        CVAdjustOrder(cv_mem, -1)
        @warn("L") = @warn("q")
        @warn("q") -= 1
        @warn("qwait") = @warn("L")
        CVRescale(cv_mem)
        return FALSE
    end
    @warn("eta") = @warn("MAX(ETAMIN,hmin/ABS(h))")
    @warn("h") *= @warn("eta")
    @warn("hscale") = @warn("h")
    @warn("qwait") = LONG_WAIT
    (@warn("N"), @warn("tn"), (@warn("zn"))[1], @warn("tempv"), @warn("f_data"))
    @warn("nfe") += 1
    N_VScale(@warn("h"), @warn("tempv"), (@warn("zn"))[2])
    return FALSE
end
function CVCompleteStep(cv_mem)
    local i
    local j
    @warn("nst") += 1
    @warn("hu") = @warn("h")
    @warn("qu") = @warn("q")
    (@warn("tau"))[i + 1] = (@warn("tau"))[(i - 1) + 1]
    if @warn("q") == 1 && @warn("nst") > 1
        (@warn("tau"))[3] = (@warn("tau"))[2]
    end
    (@warn("tau"))[2] = @warn("h")
    N_VLinearSum((@warn("l"))[j + 1], @warn("acor"), ONE, (@warn("zn"))[j + 1], (@warn("zn"))[j + 1])
    @warn("qwait") -= 1
    if @warn("qwait") == 1 && @warn("q") != @warn("qmax")
        N_VScale(ONE, @warn("acor"), (@warn("zn"))[@warn("qmax") + 1])
        @warn("saved_tq5") = (@warn("tq"))[6]
    end
end
function CVPrepareNextStep(cv_mem, dsm)
    local etaqm1
    local etaq
    local etaqp1
    if @warn("etamax") == ONE
        @warn("qwait") = @warn("MAX(qwait,2)")
        @warn("qprime") = @warn("q")
        @warn("hprime") = @warn("h")
        @warn("eta") = ONE
        @warn("etamax") = if @warn("nst") <= SMALL_NST
                ETAMX2
            else
                ETAMX3
            end
        N_VScale(ONE / (@warn("tq"))[3], @warn("acor"), @warn("acor"))
        return
    end
    etaq = ONE / (RPowerR(BIAS2 * dsm, ONE / @warn("L")) + ADDON)
    if @warn("qwait") != 0
        @warn("eta") = etaq
        @warn("qprime") = @warn("q")
        CVSetEta(cv_mem)
        return
    end
    @warn("qwait") = 2
    etaqm1 = CVComputeEtaqm1(cv_mem)
    etaqp1 = CVComputeEtaqp1(cv_mem)
    CVChooseEta(cv_mem, etaqm1, etaq, etaqp1)
    CVSetEta(cv_mem)
end
function CVSetEta(cv_mem)
    if @warn("eta") < THRESH
        @warn("eta") = ONE
        @warn("hprime") = @warn("h")
    else
        @warn("eta") = @warn("MIN(eta,etamax)")
        @warn("eta") /= @warn("MAX(ONE,ABS(h)*hmax_inv*eta)")
        @warn("hprime") = @warn("h") * @warn("eta")
    end
    @warn("etamax") = if @warn("nst") <= SMALL_NST
            ETAMX2
        else
            ETAMX3
        end
    N_VScale(ONE / (@warn("tq"))[3], @warn("acor"), @warn("acor"))
end
function CVComputeEtaqm1(cv_mem)
    local etaqm1
    local ddn
    etaqm1 = ZERO
    if @warn("q") > 1
        ddn = N_VWrmsNorm((@warn("zn"))[@warn("q") + 1], @warn("ewt")) / (@warn("tq"))[2]
        etaqm1 = ONE / (RPowerR(BIAS1 * ddn, ONE / @warn("q")) + ADDON)
    end
    return etaqm1
end
function CVComputeEtaqp1(cv_mem)
    local etaqp1
    local dup
    local cquot
    etaqp1 = ZERO
    if @warn("q") != @warn("qmax")
        cquot = ((@warn("tq"))[6] / @warn("saved_tq5")) * RPowerI(@warn("h") / (@warn("tau"))[3], @warn("L"))
        N_VLinearSum(-cquot, (@warn("zn"))[@warn("qmax") + 1], ONE, @warn("acor"), @warn("tempv"))
        dup = N_VWrmsNorm(@warn("tempv"), @warn("ewt")) / (@warn("tq"))[4]
        etaqp1 = ONE / (RPowerR(BIAS3 * dup, ONE / (@warn("L") + 1)) + ADDON)
    end
    return etaqp1
end
function CVChooseEta(cv_mem, etaqm1, etaq, etaqp1)
    local etam
    etam = @warn("MAX(etaqm1,MAX(etaq,etaqp1))")
    if etam < THRESH
        @warn("eta") = ONE
        @warn("qprime") = @warn("q")
        return
    end
    if etam == etaq
        @warn("eta") = etaq
        @warn("qprime") = @warn("q")
    elseif etam == etaqm1
        @warn("eta") = etaqm1
        @warn("qprime") = @warn("q") - 1
    else
        @warn("eta") = etaqp1
        @warn("qprime") = @warn("q") + 1
        N_VScale(ONE, @warn("acor"), (@warn("zn"))[@warn("qmax") + 1])
    end
end
function CVHandleFailure(cv_mem, kflag)
    N_VProd(@warn("acor"), @warn("ewt"), @warn("tempv"))
    N_VAbs(@warn("tempv"), @warn("tempv"))
    @switch kflag begin
            @case ?_?(REP_ERR_FAIL)
            fprintf(@warn("errfp"), MSG_ERR_FAILS, @warn("tn"), @warn("h"))
            return ERR_FAILURE
            @case ?_?(REP_CONV_FAIL)
            fprintf(@warn("errfp"), MSG_CONV_FAILS, @warn("tn"), @warn("h"))
            return CONV_FAILURE
            @case ?_?(SETUP_FAILED)
            fprintf(@warn("errfp"), MSG_SETUP_FAILED, @warn("tn"))
            return SETUP_FAILURE
            @case ?_?(SOLVE_FAILED)
            fprintf(@warn("errfp"), MSG_SOLVE_FAILED, @warn("tn"))
            return SOLVE_FAILURE
        end
end
