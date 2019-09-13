function idamax(n, dx, incx)
    local dmax
    local xmag
    local i
    local ii
    local xindex
    xindex = 0
    if n <= 0
        return xindex
    end
    xindex = 1
    if n <= 1 || incx <= 0
        return xindex
    end
    if incx != 1
        dmax = fabs(dx[1])
        ii = 2
        begin
            xmag = fabs(dx[i + 1])
            if xmag > dmax
                xindex = ii
                dmax = xmag
            end
            ii += 1
        end
        return xindex
    end
    dmax = fabs(dx[1])
    begin
        xmag = fabs(dx[i + 1])
        if xmag > dmax
            xindex = i
            dmax = xmag
        end
    end
    return xindex
end
function dscal(n, da, dx, incx)
    local m
    local i
    if n <= 0
        return
    end
    if incx != 1
        dx[i + 1] = da * dx[i + 1]
        return
    end
    m = n % 5
    if m != 0
        dx[i + 1] = da * dx[i + 1]
        if n < 5
            return
        end
    end
    begin
        dx[i + 1] = da * dx[i + 1]
        dx[(i + 1) + 1] = da * dx[(i + 1) + 1]
        dx[(i + 2) + 1] = da * dx[(i + 2) + 1]
        dx[(i + 3) + 1] = da * dx[(i + 3) + 1]
        dx[(i + 4) + 1] = da * dx[(i + 4) + 1]
    end
    return
end
function ddot(n, dx, incx, dy, incy)
    local dotprod
    local ix
    local iy
    local i
    local m
    dotprod = 0.0
    if n <= 0
        return dotprod
    end
    if incx != incy || incx < 1
        ix = 1
        iy = 1
        if incx < 0
            ix = (-n + 1) * incx + 1
        end
        if incy < 0
            iy = (-n + 1) * incy + 1
        end
        begin
            dotprod = dotprod + dx[ix + 1] * dy[iy + 1]
            ix = ix + incx
            iy = iy + incy
        end
        return dotprod
    end
    if incx == 1
        m = n % 5
        if m != 0
            dotprod = dotprod + dx[i + 1] * dy[i + 1]
            if n < 5
                return dotprod
            end
        end
        dotprod = ((((dotprod + dx[i + 1] * dy[i + 1]) + dx[(i + 1) + 1] * dy[(i + 1) + 1]) + dx[(i + 2) + 1] * dy[(i + 2) + 1]) + dx[(i + 3) + 1] * dy[(i + 3) + 1]) + dx[(i + 4) + 1] * dy[(i + 4) + 1]
        return dotprod
    end
    dotprod = dotprod + dx[i + 1] * dy[i + 1]
    return dotprod
end
function daxpy(n, da, dx, incx, dy, incy)
    local ix
    local iy
    local i
    local m
    if n < 0 || da == 0.0
        return
    end
    if incx != incy || incx < 1
        ix = 1
        iy = 1
        if incx < 0
            ix = (-n + 1) * incx + 1
        end
        if incy < 0
            iy = (-n + 1) * incy + 1
        end
        begin
            dy[iy + 1] = dy[iy + 1] + da * dx[ix + 1]
            ix = ix + incx
            iy = iy + incy
        end
        return
    end
    if incx == 1
        m = n % 4
        if m != 0
            dy[i + 1] = dy[i + 1] + da * dx[i + 1]
            if n < 4
                return
            end
        end
        begin
            dy[i + 1] = dy[i + 1] + da * dx[i + 1]
            dy[(i + 1) + 1] = dy[(i + 1) + 1] + da * dx[(i + 1) + 1]
            dy[(i + 2) + 1] = dy[(i + 2) + 1] + da * dx[(i + 2) + 1]
            dy[(i + 3) + 1] = dy[(i + 3) + 1] + da * dx[(i + 3) + 1]
        end
        return
    end
    dy[i + 1] = da * dx[i + 1] + dy[i + 1]
    return
end
function dgesl(a, n, ipvt, b, job)
    local nm1
    local k
    local j
    local t
    nm1 = n - 1
    if job == 0
        begin
            t = ddot(k - 1, a[k + 1], 1, b, 1)
            b[k + 1] = (b[k + 1] - t) / (a[k + 1])[k + 1]
        end
        begin
            b[k + 1] = b[k + 1] + ddot(n - k, a[k + 1] + k, 1, b + k, 1)
            j = ipvt[k + 1]
            if j != k
                t = b[j + 1]
                b[j + 1] = b[k + 1]
                b[k + 1] = t
            end
        end
        return
    end
    begin
        j = ipvt[k + 1]
        t = b[j + 1]
        if j != k
            b[j + 1] = b[k + 1]
            b[k + 1] = t
        end
        daxpy(n - k, t, a[k + 1] + k, 1, b + k, 1)
    end
    begin
        b[k + 1] = b[k + 1] / (a[k + 1])[k + 1]
        t = -(b[k + 1])
        daxpy(k - 1, t, a[k + 1], 1, b, 1)
    end
end
function dgefa(a, n, ipvt, info)
    local j
    local k
    local i
    local t
    info[] = 0
    begin
        j = (idamax((n - k) + 1, (a[k + 1] + k) - 1, 1) + k) - 1
        ipvt[k + 1] = j
        if (a[k + 1])[j + 1] == 0.0
            info[] = k
            continue
        end
        if j != k
            t = (a[k + 1])[j + 1]
            (a[k + 1])[j + 1] = (a[k + 1])[k + 1]
            (a[k + 1])[k + 1] = t
        end
        t = -1.0 / (a[k + 1])[k + 1]
        dscal(n - k, t, a[k + 1] + k, 1)
        begin
            t = (a[i + 1])[j + 1]
            if j != k
                (a[i + 1])[j + 1] = (a[i + 1])[k + 1]
                (a[i + 1])[k + 1] = t
            end
            daxpy(n - k, t, a[k + 1] + k, 1, a[i + 1] + k, 1)
        end
    end
    ipvt[n + 1] = n
    if (a[n + 1])[n + 1] == 0.0
        info[] = n
    end
end
function stoda(neq, y, f, _data)
    local corflag
    local orderflag
    local i
    local i1
    local j
    local m
    local ncf
    local del
    local delp
    local dsm
    local dup
    local exup
    local r
    local rh
    local rhup
    local told
    local pdh
    local pnorm
    kflag = 0
    told = tn
    ncf = 0
    ierpj = 0
    iersl = 0
    jcur = 0
    delp = 0.0
    if jstart == 0
        lmax = maxord + 1
        nq = 1
        l = 2
        ialth = 2
        rmax = 10000.0
        rc = 0.0
        el0 = 1.0
        crate = 0.7
        hold = h
        nslp = 0
        ipup = miter
        icount = 20
        irflag = 0
        pdest = 0.0
        pdlast = 0.0
        ratio = 5.0
        cfode(2)
        cm2[i + 1] = (tesco[i + 1])[2] * (elco[i + 1])[(i + 1) + 1]
        cfode(1)
        cm1[i + 1] = (tesco[i + 1])[2] * (elco[i + 1])[(i + 1) + 1]
        resetcoeff()
    end
    if jstart == -1
        ipup = miter
        lmax = maxord + 1
        if ialth == 1
            ialth = 2
        end
        if meth != mused
            cfode(meth)
            ialth = l
            resetcoeff()
        end
        if h != hold
            rh = h / hold
            h = hold
            @c scaleh(&rh, &pdh)
        end
    end
    if jstart == -2
        if h != hold
            rh = h / hold
            h = hold
            @c scaleh(&rh, &pdh)
        end
    end
    while 1
        while 1
            if fabs(rc - 1.0) > ccmax
                ipup = miter
            end
            if nst >= nslp + msbp
                ipup = miter
            end
            tn += h
            begin
                yp1 = yh[i1 + 1]
                yp2 = yh[(i1 + 1) + 1]
                yp1[i + 1] += yp2[i + 1]
            end
            pnorm = vmnorm(n, yh[1], ewt)
            @c correction(neq, y, f, &corflag, pnorm, &del, &delp, &told, &ncf, &rh, &m, _data)
            if corflag == 0
                break
            end
            if corflag == 1
                rh = @warn("max(rh,hmin/fabs(h))")
                @c scaleh(&rh, &pdh)
                continue
            end
            if corflag == 2
                kflag = -2
                hold = h
                jstart = 1
                return
            end
        end
        jcur = 0
        if m == 0
            dsm = del / (tesco[nq + 1])[2]
        end
        if m > 0
            dsm = vmnorm(n, acor, ewt) / (tesco[nq + 1])[2]
        end
        if dsm <= 1.0
            kflag = 0
            nst += 1
            hu = h
            nqu = nq
            mused = meth
            begin
                yp1 = yh[j + 1]
                r = el[j + 1]
                yp1[i + 1] += r * acor[i + 1]
            end
            icount -= 1
            if icount < 0
                @c methodswitch(dsm, pnorm, &pdh, &rh)
                if meth != mused
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                    rmax = 10.0
                    endstoda()
                    break
                end
            end
            ialth -= 1
            if ialth == 0
                rhup = 0.0
                if l != lmax
                    yp1 = yh[lmax + 1]
                    savf[i + 1] = acor[i + 1] - yp1[i + 1]
                    dup = vmnorm(n, savf, ewt) / (tesco[nq + 1])[3]
                    exup = 1.0 / Cdouble(l + 1)
                    rhup = 1.0 / (1.4 * pow(dup, exup) + 1.4e-6)
                end
                @c orderswitch(&rhup, dsm, &pdh, &rh, &orderflag)
                if orderflag == 0
                    endstoda()
                    break
                end
                if orderflag == 1
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                    rmax = 10.0
                    endstoda()
                    break
                end
                if orderflag == 2
                    resetcoeff()
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                    rmax = 10.0
                    endstoda()
                    break
                end
            end
            if ialth > 1 || l == lmax
                endstoda()
                break
            end
            yp1 = yh[lmax + 1]
            yp1[i + 1] = acor[i + 1]
            endstoda()
            break
        else
            kflag -= 1
            tn = told
            begin
                yp1 = yh[i1 + 1]
                yp2 = yh[(i1 + 1) + 1]
                yp1[i + 1] -= yp2[i + 1]
            end
            rmax = 2.0
            if fabs(h) <= hmin * 1.00001
                kflag = -1
                hold = h
                jstart = 1
                break
            end
            if kflag > -3
                rhup = 0.0
                @c orderswitch(&rhup, dsm, &pdh, &rh, &orderflag)
                if orderflag == 1 || orderflag == 0
                    if orderflag == 0
                        rh = @warn("min(rh,0.2)")
                    end
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                end
                if orderflag == 2
                    resetcoeff()
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                end
                continue
            else
                if kflag == -10
                    kflag = -1
                    hold = h
                    jstart = 1
                    break
                else
                    rh = 0.1
                    rh = @warn("max(hmin/fabs(h),rh)")
                    h *= rh
                    yp1 = yh[1]
                    y[i + 1] = yp1[i + 1]
                    (tn, y + 1, savf + 1, _data)
                    nfe += 1
                    yp1 = yh[2]
                    yp1[i + 1] = h * savf[i + 1]
                    ipup = miter
                    ialth = 5
                    if nq == 1
                        continue
                    end
                    nq = 1
                    l = 2
                    resetcoeff()
                    continue
                end
            end
        end
    end
end
function correction(neq, y, f, corflag, pnorm, del, delp, told, ncf, rh, m, _data)
    local i
    local rm
    local rate
    local dcon
    m[] = 0
    corflag[] = 0
    rate = 0.0
    del[] = 0.0
    yp1 = yh[1]
    y[i + 1] = yp1[i + 1]
    (tn, y + 1, savf + 1, _data)
    nfe += 1
    while 1
        if m[] == 0
            if ipup > 0
                prja(neq, y, f, _data)
                ipup = 0
                rc = 1.0
                nslp = nst
                crate = 0.7
                if ierpj != 0
                    corfailure(told, rh, ncf, corflag)
                    return
                end
            end
            acor[i + 1] = 0.0
        end
        if miter == 0
            yp1 = yh[2]
            begin
                savf[i + 1] = h * savf[i + 1] - yp1[i + 1]
                y[i + 1] = savf[i + 1] - acor[i + 1]
            end
            del[] = vmnorm(n, y, ewt)
            yp1 = yh[1]
            begin
                y[i + 1] = yp1[i + 1] + el[1] * savf[i + 1]
                acor[i + 1] = savf[i + 1]
            end
        else
            yp1 = yh[2]
            y[i + 1] = h * savf[i + 1] - (yp1[i + 1] + acor[i + 1])
            solsy(y)
            del[] = vmnorm(n, y, ewt)
            yp1 = yh[1]
            begin
                acor[i + 1] += y[i + 1]
                y[i + 1] = yp1[i + 1] + el[1] * acor[i + 1]
            end
        end
        if del[] <= (100.0pnorm) * ETA
            break
        end
        if m[] != 0 || meth != 1
            if m[] != 0
                rm = 1024.0
                if del[] <= 1024.0 * delp[]
                    rm = del[] / delp[]
                end
                rate = @warn("max(rate,rm)")
                crate = @warn("max(0.2*crate,rm)")
            end
            dcon = (del[] * @warn("min(1.,1.5*crate)")) / ((tesco[nq + 1])[2] * conit)
            if dcon <= 1.0
                pdest = @warn("max(pdest,rate/fabs(h*el[1]))")
                if pdest != 0.0
                    pdlast = pdest
                end
                break
            end
        end
        (m[])[]
        if m[] == maxcor || m[] >= 2 && del[] > 2.0 * delp[]
            if miter == 0 || jcur == 1
                corfailure(told, rh, ncf, corflag)
                return
            end
            ipup = miter
            m[] = 0
            rate = 0.0
            del[] = 0.0
            yp1 = yh[1]
            y[i + 1] = yp1[i + 1]
            (tn, y + 1, savf + 1, _data)
            nfe += 1
        else
            delp[] = del[]
            (tn, y + 1, savf + 1, _data)
            nfe += 1
        end
    end
end
function prja(neq, y, f, _data)
    local i
    local ier
    local j
    local fac
    local hl0
    local r
    local r0
    local yj
    nje += 1
    ierpj = 0
    jcur = 1
    hl0 = h * el0
    if miter != 2
        fprintf(__stderrp, "[prja] miter != 2\n")
        return
    end
    if miter == 2
        fac = vmnorm(n, savf, ewt)
        r0 = (((1000.0 * fabs(h)) * ETA) * Cdouble(n)) * fac
        if r0 == 0.0
            r0 = 1.0
        end
        begin
            yj = y[j + 1]
            r = @warn("max(sqrteta*fabs(yj),r0/ewt[j])")
            y[j + 1] += r
            fac = -hl0 / r
            (tn, y + 1, acor + 1, _data)
            (wm[i + 1])[j + 1] = (acor[i + 1] - savf[i + 1]) * fac
            y[j + 1] = yj
        end
        nfe += n
        pdnorm = fnorm(n, wm, ewt) / fabs(hl0)
        (wm[i + 1])[i + 1] += 1.0
        @c dgefa(wm, n, ipvt, &ier)
        if ier != 0
            ierpj = 1
        end
        return
    end
end
function terminate(istate)
    if illin == 5
        fprintf(__stderrp, "[lsoda] repeated occurrence of illegal input. run aborted.. apparent infinite loop\n")
    else
        illin += 1
        istate[] = -3
    end
end
function terminate2(y, t)
    local i
    yp1 = yh[1]
    y[i + 1] = yp1[i + 1]
    t[] = tn
    illin = 0
    freevectors()
    return
end
function successreturn(y, t, itask, ihit, tcrit, istate)
    local i
    yp1 = yh[1]
    y[i + 1] = yp1[i + 1]
    t[] = tn
    if itask == 4 || itask == 5
        elseif ihit
            t[] = tcrit
        end
    end
    istate[] = 2
    illin = 0
    freevectors()
end
function freevectors()
end
function _freevectors()
    local i
    if wm
        free(wm[i + 1])
    end
    if yh
        free(yh[i + 1])
    end
    free(yh)
    free(wm)
    free(ewt)
    free(savf)
    free(acor)
    free(ipvt)
    g_nyh = (g_lenyh = 0)
    yh = 0
    wm = 0
    ewt = 0
    savf = 0
    acor = 0
    ipvt = 0
end
function ewset(itol, rtol, atol, ycur)
    local i
    @switch itol begin
            @case 1
            ewt[i + 1] = rtol[1] * fabs(ycur[i + 1]) + atol[1]
            break
            @case 2
            ewt[i + 1] = rtol[1] * fabs(ycur[i + 1]) + atol[i + 1]
            break
            @case 3
            ewt[i + 1] = rtol[i + 1] * fabs(ycur[i + 1]) + atol[1]
            break
            @case 4
            ewt[i + 1] = rtol[i + 1] * fabs(ycur[i + 1]) + atol[i + 1]
            break
        end
end
function resetcoeff()
    local i
    local ep1
    ep1 = elco[nq + 1]
    el[i + 1] = ep1[i + 1]
    rc = (rc * el[1]) / el0
    el0 = el[1]
    conit = 0.5 / Cdouble(nq + 2)
end
function solsy(y)
    iersl = 0
    if miter != 2
        printf("solsy -- miter != 2\n")
        return
    end
    if miter == 2
        dgesl(wm, n, ipvt, y, 0)
    end
    return
end
function endstoda()
    local r
    local i
    r = 1.0 / (tesco[nqu + 1])[2]
    acor[i + 1] *= r
    hold = h
    jstart = 1
end
function orderswitch(rhup, dsm, pdh, rh, orderflag)
    local newq
    local i
    local exsm
    local rhdn
    local rhsm
    local ddn
    local exdn
    local r
    orderflag[] = 0
    exsm = 1.0 / Cdouble(l)
    rhsm = 1.0 / (1.2 * pow(dsm, exsm) + 1.2e-6)
    rhdn = 0.0
    if nq != 1
        ddn = vmnorm(n, yh[l + 1], ewt) / (tesco[nq + 1])[1]
        exdn = 1.0 / Cdouble(nq)
        rhdn = 1.0 / (1.3 * pow(ddn, exdn) + 1.3e-6)
    end
    if meth == 1
        pdh[] = @warn("max(fabs(h)*pdlast,0.000001)")
        if l < lmax
            rhup[] = @warn("min(*rhup,sm1[l]/*pdh)")
        end
        rhsm = @warn("min(rhsm,sm1[nq]/*pdh)")
        if nq > 1
            rhdn = @warn("min(rhdn,sm1[nq-1]/*pdh)")
        end
        pdest = 0.0
    end
    if rhsm >= rhup[]
        if rhsm >= rhdn
            newq = nq
            rh[] = rhsm
        else
            newq = nq - 1
            rh[] = rhdn
            if kflag < 0 && rh[] > 1.0
                rh[] = 1.0
            end
        end
    else
        if rhup[] <= rhdn
            newq = nq - 1
            rh[] = rhdn
            if kflag < 0 && rh[] > 1.0
                rh[] = 1.0
            end
        else
            rh[] = rhup[]
            if rh[] >= 1.1
                r = el[l + 1] / Cdouble(l)
                nq = l
                l = nq + 1
                yp1 = yh[l + 1]
                yp1[i + 1] = acor[i + 1] * r
                orderflag[] = 2
                return
            else
                ialth = 3
                return
            end
        end
    end
    if meth == 1
        if (rh[] * pdh[]) * 1.00001 < sm1[newq + 1]
            elseif kflag == 0 && rh[] < 1.1
                ialth = 3
                return
            end
        end
    else
        if kflag == 0 && rh[] < 1.1
            ialth = 3
            return
        end
    end
    if kflag <= -2
        rh[] = @warn("min(*rh,0.2)")
    end
    if newq == nq
        orderflag[] = 1
        return
    end
    nq = newq
    l = nq + 1
    orderflag[] = 2
end
function intdy(t, k, dky, iflag)
    local i
    local ic
    local j
    local jj
    local jp1
    local c
    local r
    local s
    local tp
    iflag[] = 0
    if k < 0 || k > nq
        fprintf(__stderrp, "[intdy] k = %d illegal\n", k)
        iflag[] = -1
        return
    end
    tp = (tn - hu) - (100.0ETA) * (tn + hu)
    if (t - tp) * (t - tn) > 0.0
        fprintf(__stderrp, "intdy -- t = %g illegal. t not in interval tcur - hu to tcur\n", t)
        iflag[] = -2
        return
    end
    s = (t - tn) / h
    ic = 1
    ic *= jj
    c = Cdouble(ic)
    yp1 = yh[l + 1]
    dky[i + 1] = c * yp1[i + 1]
    begin
        jp1 = j + 1
        ic = 1
        ic *= jj
        c = Cdouble(ic)
        yp1 = yh[jp1 + 1]
        dky[i + 1] = c * yp1[i + 1] + s * dky[i + 1]
    end
    if k == 0
        return
    end
    r = pow(h, Cdouble(-k))
    dky[i + 1] *= r
end
function corfailure(told, rh, ncf, corflag)
    local j
    local i1
    local i
    ncf += 1
    rmax = 2.0
    tn = told[]
    begin
        yp1 = yh[i1 + 1]
        yp2 = yh[(i1 + 1) + 1]
        yp1[i + 1] -= yp2[i + 1]
    end
    if fabs(h) <= hmin * 1.00001 || ncf[] == mxncf
        corflag[] = 2
        return
    end
    corflag[] = 1
    rh[] = 0.25
    ipup = miter
end
function methodswitch(dsm, pnorm, pdh, rh)
    local lm1
    local lm1p1
    local lm2
    local lm2p1
    local nqm1
    local nqm2
    local rh1
    local rh2
    local rh1it
    local exm2
    local dm2
    local exm1
    local dm1
    local alpha
    local exsm
    if meth == 1
        if nq > 5
            return
        end
        if dsm <= (100.0pnorm) * ETA || pdest == 0.0
            if irflag == 0
                return
            end
            rh2 = 2.0
            nqm2 = @warn("min(nq,mxords)")
        else
            exsm = 1.0 / Cdouble(l)
            rh1 = 1.0 / (1.2 * pow(dsm, exsm) + 1.2e-6)
            rh1it = 2.0rh1
            pdh[] = pdlast * fabs(h)
            if pdh[] * rh1 > 1.0e-5
                rh1it = sm1[nq + 1] / pdh[]
            end
            rh1 = @warn("min(rh1,rh1it)")
            if nq > mxords
                nqm2 = mxords
                lm2 = mxords + 1
                exm2 = 1.0 / Cdouble(lm2)
                lm2p1 = lm2 + 1
                dm2 = vmnorm(n, yh[lm2p1 + 1], ewt) / cm2[mxords + 1]
                rh2 = 1.0 / (1.2 * pow(dm2, exm2) + 1.2e-6)
            else
                dm2 = dsm * (cm1[nq + 1] / cm2[nq + 1])
                rh2 = 1.0 / (1.2 * pow(dm2, exsm) + 1.2e-6)
                nqm2 = nq
            end
            if rh2 < ratio * rh1
                return
            end
        end
        rh[] = rh2
        icount = 20
        meth = 2
        miter = jtyp
        pdlast = 0.0
        nq = nqm2
        l = nq + 1
        return
    end
    exsm = 1.0 / Cdouble(l)
    if mxordn < nq
        nqm1 = mxordn
        lm1 = mxordn + 1
        exm1 = 1.0 / Cdouble(lm1)
        lm1p1 = lm1 + 1
        dm1 = vmnorm(n, yh[lm1p1 + 1], ewt) / cm1[mxordn + 1]
        rh1 = 1.0 / (1.2 * pow(dm1, exm1) + 1.2e-6)
    else
        dm1 = dsm * (cm2[nq + 1] / cm1[nq + 1])
        rh1 = 1.0 / (1.2 * pow(dm1, exsm) + 1.2e-6)
        nqm1 = nq
        exm1 = exsm
    end
    rh1it = 2.0rh1
    pdh[] = pdnorm * fabs(h)
    if pdh[] * rh1 > 1.0e-5
        rh1it = sm1[nqm1 + 1] / pdh[]
    end
    rh1 = @warn("min(rh1,rh1it)")
    rh2 = 1.0 / (1.2 * pow(dsm, exsm) + 1.2e-6)
    if rh1 * ratio < 5.0rh2
        return
    end
    alpha = @warn("max(0.001,rh1)")
    dm1 *= pow(alpha, exm1)
    if dm1 <= (1000.0ETA) * pnorm
        return
    end
    rh[] = rh1
    icount = 20
    meth = 1
    miter = 0
    pdlast = 0.0
    nq = nqm1
    l = nq + 1
end
function cfode(meth)
    local i
    local nq
    local nqm1
    local nqp1
    local agamq
    local fnq
    local fnqm1
    pc = 13
    local pint
    local ragq
    local rqfac
    local rq1fac
    local tsign
    local xpin
    if meth == 1
        (elco[1])[1] = 1.0
        (elco[1])[2] = 1.0
        (tesco[1])[1] = 0.0
        (tesco[1])[2] = 2.0
        (tesco[2])[1] = 1.0
        (tesco[12])[3] = 0.0
        pc[1] = 1.0
        rqfac = 1.0
        begin
            rq1fac = rqfac
            rqfac = rqfac / Cdouble(nq)
            nqm1 = nq - 1
            fnqm1 = Cdouble(nqm1)
            nqp1 = nq + 1
            pc[nq + 1] = 0.0
            pc[i + 1] = pc[(i - 1) + 1] + fnqm1 * pc[i + 1]
            pc[1] = fnqm1 * pc[1]
            pint = pc[1]
            xpin = pc[1] / 2.0
            tsign = 1.0
            begin
                tsign = -tsign
                pint += (tsign * pc[i + 1]) / Cdouble(i)
                xpin += (tsign * pc[i + 1]) / Cdouble(i + 1)
            end
            (elco[nq + 1])[1] = pint * rq1fac
            (elco[nq + 1])[2] = 1.0
            (elco[nq + 1])[(i + 1) + 1] = (rq1fac * pc[i + 1]) / Cdouble(i)
            agamq = rqfac * xpin
            ragq = 1.0 / agamq
            (tesco[nq + 1])[2] = ragq
            if nq < 12
                (tesco[nqp1 + 1])[1] = (ragq * rqfac) / Cdouble(nqp1)
            end
            (tesco[nqm1 + 1])[3] = ragq
        end
        return
    end
    pc[1] = 1.0
    rq1fac = 1.0
    begin
        fnq = Cdouble(nq)
        nqp1 = nq + 1
        pc[nqp1 + 1] = 0.0
        pc[i + 1] = pc[(i - 1) + 1] + fnq * pc[i + 1]
        pc[1] *= fnq
        (elco[nq + 1])[i + 1] = pc[i + 1] / pc[2]
        (elco[nq + 1])[2] = 1.0
        (tesco[nq + 1])[1] = rq1fac
        (tesco[nq + 1])[2] = Cdouble(nqp1) / (elco[nq + 1])[1]
        (tesco[nq + 1])[3] = Cdouble(nq + 2) / (elco[nq + 1])[1]
        rq1fac /= fnq
    end
    return
end
function scaleh(rh, pdh)
    local r
    local j
    local i
    rh[] = @warn("min(*rh,rmax)")
    rh[] = rh[] / @warn("max(1.,fabs(h)*hmxi**rh)")
    if meth == 1
        irflag = 0
        pdh[] = @warn("max(fabs(h)*pdlast,0.000001)")
        if (rh[] * pdh[]) * 1.00001 >= sm1[nq + 1]
            rh[] = sm1[nq + 1] / pdh[]
            irflag = 1
        end
    end
    r = 1.0
    begin
        r *= rh[]
        yp1 = yh[j + 1]
        yp1[i + 1] *= r
    end
    h *= rh[]
    rc *= rh[]
    ialth = l
end
function fnorm(n, a, w)
    local i
    local j
    local an
    local sum
    local ap1
    an = 0.0
    begin
        sum = 0.0
        ap1 = a[i + 1]
        sum += fabs(ap1[j + 1]) / w[j + 1]
        an = @warn("max(an,sum*w[i])")
    end
    return an
end
function vmnorm(n, v, w)
    local i
    local vm
    vm = 0.0
    vm = @warn("max(vm,fabs(v[i])*w[i])")
    return vm
end
function terminate(istate)
    if illin == 5
        fprintf(__stderrp, "[lsoda] repeated occurrence of illegal input. run aborted.. apparent infinite loop\n")
    else
        illin += 1
        istate[] = -3
    end
end
function terminate2(y, t)
    local i
    yp1 = yh[1]
    y[i + 1] = yp1[i + 1]
    t[] = tn
    illin = 0
    freevectors()
    return
end
function successreturn(y, t, itask, ihit, tcrit, istate)
    local i
    yp1 = yh[1]
    y[i + 1] = yp1[i + 1]
    t[] = tn
    if itask == 4 || itask == 5
        elseif ihit
            t[] = tcrit
        end
    end
    istate[] = 2
    illin = 0
    freevectors()
end
function lsoda(f, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, jt, iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9, rwork1, rwork5, rwork6, rwork7, _data)
    mxstp0 = 500
    mxhnl0 = 10
    local i
    local iflag
    local lenyh
    local ihit
    local atoli
    local ayi
    local big
    local h0
    local hmax
    local hmx
    local rh
    local rtoli
    local tcrit
    local tdist
    local tnext
    local tol
    local tolsf
    local tp
    local size
    local sum
    local w0
    if istate[] == 1
        _freevectors()
    end
    if istate[] < 1 || istate[] > 3
        fprintf(__stderrp, "[lsoda] illegal istate = %d\n", istate[])
        terminate(istate)
        return
    end
    if itask < 1 || itask > 5
        fprintf(__stderrp, "[lsoda] illegal itask = %d\n", itask)
        terminate(istate)
        return
    end
    if init == 0 && (istate[] == 2 || istate[] == 3)
        fprintf(__stderrp, "[lsoda] istate > 1 but lsoda not initialized\n")
        terminate(istate)
        return
    end
    if istate[] == 1
        init = 0
        if tout == t[]
            ntrep += 1
            if ntrep < 5
                return
            end
            fprintf(__stderrp, "[lsoda] repeated calls with istate = 1 and tout = t. run aborted.. apparent infinite loop\n")
            return
        end
    end
    if istate[] == 1 || istate[] == 3
        ntrep = 0
        if neq <= 0
            fprintf(__stderrp, "[lsoda] neq = %d is less than 1\n", neq)
            terminate(istate)
            return
        end
        if istate[] == 3 && neq > n
            fprintf(__stderrp, "[lsoda] istate = 3 and neq increased\n")
            terminate(istate)
            return
        end
        n = neq
        if itol < 1 || itol > 4
            fprintf(__stderrp, "[lsoda] itol = %d illegal\n", itol)
            terminate(istate)
            return
        end
        if iopt < 0 || iopt > 1
            fprintf(__stderrp, "[lsoda] iopt = %d illegal\n", iopt)
            terminate(istate)
            return
        end
        if (jt == 3 || jt < 1) || jt > 5
            fprintf(__stderrp, "[lsoda] jt = %d illegal\n", jt)
            terminate(istate)
            return
        end
        jtyp = jt
        if jt > 2
            ml = iwork1
            mu = iwork2
            if ml < 0 || ml >= n
                fprintf(__stderrp, "[lsoda] ml = %d not between 1 and neq\n", ml)
                terminate(istate)
                return
            end
            if mu < 0 || mu >= n
                fprintf(__stderrp, "[lsoda] mu = %d not between 1 and neq\n", mu)
                terminate(istate)
                return
            end
        end
        if iopt == 0
            ixpr = 0
            mxstep = mxstp0
            mxhnil = mxhnl0
            hmxi = 0.0
            hmin = 0.0
            if istate[] == 1
                h0 = 0.0
                mxordn = mord[1]
                mxords = mord[2]
            end
        else
            ixpr = iwork5
            if ixpr < 0 || ixpr > 1
                fprintf(__stderrp, "[lsoda] ixpr = %d is illegal\n", ixpr)
                terminate(istate)
                return
            end
            mxstep = iwork6
            if mxstep < 0
                fprintf(__stderrp, "[lsoda] mxstep < 0\n")
                terminate(istate)
                return
            end
            if mxstep == 0
                mxstep = mxstp0
            end
            mxhnil = iwork7
            if mxhnil < 0
                fprintf(__stderrp, "[lsoda] mxhnil < 0\n")
                terminate(istate)
                return
            end
            if istate[] == 1
                h0 = rwork5
                mxordn = iwork8
                if mxordn < 0
                    fprintf(__stderrp, "[lsoda] mxordn = %d is less than 0\n", mxordn)
                    terminate(istate)
                    return
                end
                if mxordn == 0
                    mxordn = 100
                end
                mxordn = @warn("min(mxordn,mord[1])")
                mxords = iwork9
                if mxords < 0
                    fprintf(__stderrp, "[lsoda] mxords = %d is less than 0\n", mxords)
                    terminate(istate)
                    return
                end
                if mxords == 0
                    mxords = 100
                end
                mxords = @warn("min(mxords,mord[2])")
                if (tout - t[]) * h0 < 0.0
                    fprintf(__stderrp, "[lsoda] tout = %g behind t = %g. integration direction is given by %g\n", tout, t[], h0)
                    terminate(istate)
                    return
                end
            end
            hmax = rwork6
            if hmax < 0.0
                fprintf(__stderrp, "[lsoda] hmax < 0.\n")
                terminate(istate)
                return
            end
            hmxi = 0.0
            if hmax > 0
                hmxi = 1.0 / hmax
            end
            hmin = rwork7
            if hmin < 0.0
                fprintf(__stderrp, "[lsoda] hmin < 0.\n")
                terminate(istate)
                return
            end
        end
    end
    if istate[] == 1
        sqrteta = sqrt(ETA)
        meth = 1
        g_nyh = (nyh = n)
        g_lenyh = (lenyh = 1 + @warn("max(mxordn,mxords)"))
        yh = Ptr{Ptr{Cdouble}}(calloc(1 + lenyh, yh[]))
        if yh == @warn("NULL")
            printf("lsoda -- insufficient memory for your problem\n")
            terminate(istate)
            return
        end
        yh[i + 1] = Ptr{Cdouble}(calloc(1 + nyh, nothing))
        wm = Ptr{Ptr{Cdouble}}(calloc(1 + nyh, wm[]))
        if wm == @warn("NULL")
            free(yh)
            printf("lsoda -- insufficient memory for your problem\n")
            terminate(istate)
            return
        end
        wm[i + 1] = Ptr{Cdouble}(calloc(1 + nyh, nothing))
        ewt = Ptr{Cdouble}(calloc(1 + nyh, nothing))
        if ewt == @warn("NULL")
            free(yh)
            free(wm)
            printf("lsoda -- insufficient memory for your problem\n")
            terminate(istate)
            return
        end
        savf = Ptr{Cdouble}(calloc(1 + nyh, nothing))
        if savf == @warn("NULL")
            free(yh)
            free(wm)
            free(ewt)
            printf("lsoda -- insufficient memory for your problem\n")
            terminate(istate)
            return
        end
        acor = Ptr{Cdouble}(calloc(1 + nyh, nothing))
        if acor == @warn("NULL")
            free(yh)
            free(wm)
            free(ewt)
            free(savf)
            printf("lsoda -- insufficient memory for your problem\n")
            terminate(istate)
            return
        end
        ipvt = Ptr{Cint}(calloc(1 + nyh, nothing))
        if ipvt == @warn("NULL")
            free(yh)
            free(wm)
            free(ewt)
            free(savf)
            free(acor)
            printf("lsoda -- insufficient memory for your problem\n")
            terminate(istate)
            return
        end
    end
    if istate[] == 1 || istate[] == 3
        rtoli = rtol[1]
        atoli = atol[1]
        begin
            if itol >= 3
                rtoli = rtol[i + 1]
            end
            if itol == 2 || itol == 4
                atoli = atol[i + 1]
            end
            if rtoli < 0.0
                fprintf(__stderrp, "[lsoda] rtol = %g is less than 0.\n", rtoli)
                terminate(istate)
                freevectors()
                return
            end
            if atoli < 0.0
                fprintf(__stderrp, "[lsoda] atol = %g is less than 0.\n", atoli)
                terminate(istate)
                freevectors()
                return
            end
        end
    end
    if istate[] == 3
        jstart = -1
    end
    if istate[] == 1
        tn = t[]
        tsw = t[]
        maxord = mxordn
        if itask == 4 || itask == 5
            tcrit = rwork1
            if (tcrit - tout) * (tout - t[]) < 0.0
                fprintf(__stderrp, "[lsoda] itask = 4 or 5 and tcrit behind tout\n")
                terminate(istate)
                freevectors()
                return
            end
            if h0 != 0.0 && ((t[] + h0) - tcrit) * h0 > 0.0
                h0 = tcrit - t[]
            end
        end
        jstart = 0
        nhnil = 0
        nst = 0
        nje = 0
        nslast = 0
        hu = 0.0
        nqu = 0
        mused = 0
        miter = 0
        ccmax = 0.3
        maxcor = 3
        msbp = 20
        mxncf = 10
        (t[], y + 1, yh[2] + 1, _data)
        nfe = 1
        yp1 = yh[1]
        yp1[i + 1] = y[i + 1]
        nq = 1
        h = 1.0
        ewset(itol, rtol, atol, y)
        begin
            if ewt[i + 1] <= 0.0
                fprintf(__stderrp, "[lsoda] ewt[%d] = %g <= 0.\n", i, ewt[i + 1])
                terminate2(y, t)
                return
            end
            ewt[i + 1] = 1.0 / ewt[i + 1]
        end
        if h0 == 0.0
            tdist = fabs(tout - t[])
            w0 = @warn("max(fabs(*t),fabs(tout))")
            if tdist < (2.0ETA) * w0
                fprintf(__stderrp, "[lsoda] tout too close to t to start integration\n ")
                terminate(istate)
                freevectors()
                return
            end
            tol = rtol[1]
            if itol > 2
                tol = @warn("max(tol,rtol[i])")
            end
            if tol <= 0.0
                atoli = atol[1]
                begin
                    if itol == 2 || itol == 4
                        atoli = atol[i + 1]
                    end
                    ayi = fabs(y[i + 1])
                    if ayi != 0.0
                        tol = @warn("max(tol,atoli/ayi)")
                    end
                end
            end
            tol = @warn("max(tol,100.*ETA)")
            tol = @warn("min(tol,0.001)")
            sum = vmnorm(n, yh[2], ewt)
            sum = 1.0 / ((tol * w0) * w0) + (tol * sum) * sum
            h0 = 1.0 / sqrt(sum)
            h0 = @warn("min(h0,tdist)")
            h0 = h0 * if tout - t[] >= 0.0
                        1.0
                    else
                        -1.0
                    end
        end
        rh = fabs(h0) * hmxi
        if rh > 1.0
            h0 /= rh
        end
        h = h0
        yp1 = yh[2]
        yp1[i + 1] *= h0
    end
    if istate[] == 2 || istate[] == 3
        nslast = nst
        @switch itask begin
                @case 1
                if (tn - tout) * h >= 0.0
                    @c intdy(tout, 0, y, &iflag)
                    if iflag != 0
                        fprintf(__stderrp, "[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask, tout)
                        terminate(istate)
                        freevectors()
                        return
                    end
                    t[] = tout
                    istate[] = 2
                    illin = 0
                    freevectors()
                    return
                end
                break
                @case 2
                break
                @case 3
                tp = tn - hu * (1.0 + 100.0ETA)
                if (tp - tout) * h > 0.0
                    fprintf(__stderrp, "[lsoda] itask = %d and tout behind tcur - hu\n", itask)
                    terminate(istate)
                    freevectors()
                    return
                end
                if (tn - tout) * h < 0.0
                    break
                end
                successreturn(y, t, itask, ihit, tcrit, istate)
                return
                @case 4
                tcrit = rwork1
                if (tn - tcrit) * h > 0.0
                    fprintf(__stderrp, "[lsoda] itask = 4 or 5 and tcrit behind tcur\n")
                    terminate(istate)
                    freevectors()
                    return
                end
                if (tcrit - tout) * h < 0.0
                    fprintf(__stderrp, "[lsoda] itask = 4 or 5 and tcrit behind tout\n")
                    terminate(istate)
                    freevectors()
                    return
                end
                if (tn - tout) * h >= 0.0
                    @c intdy(tout, 0, y, &iflag)
                    if iflag != 0
                        fprintf(__stderrp, "[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask, tout)
                        terminate(istate)
                        freevectors()
                        return
                    end
                    t[] = tout
                    istate[] = 2
                    illin = 0
                    freevectors()
                    return
                end
                @case 5
                if itask == 5
                    tcrit = rwork1
                    if (tn - tcrit) * h > 0.0
                        fprintf(__stderrp, "[lsoda] itask = 4 or 5 and tcrit behind tcur\n")
                        terminate(istate)
                        freevectors()
                        return
                    end
                end
                hmx = fabs(tn) + fabs(h)
                ihit = fabs(tn - tcrit) <= (100.0ETA) * hmx
                if ihit
                    t[] = tcrit
                    successreturn(y, t, itask, ihit, tcrit, istate)
                    return
                end
                tnext = tn + h * (1.0 + 4.0ETA)
                if (tnext - tcrit) * h <= 0.0
                    break
                end
                h = (tcrit - tn) * (1.0 - 4.0ETA)
                if istate[] == 2
                    jstart = -2
                end
                break
            end
    end
    while 1
        if istate[] != 1 || nst != 0
            if nst - nslast >= mxstep
                fprintf(__stderrp, "[lsoda] %d steps taken before reaching tout\n", mxstep)
                istate[] = -1
                terminate2(y, t)
                return
            end
            ewset(itol, rtol, atol, yh[1])
            begin
                if ewt[i + 1] <= 0.0
                    fprintf(__stderrp, "[lsoda] ewt[%d] = %g <= 0.\n", i, ewt[i + 1])
                    istate[] = -6
                    terminate2(y, t)
                    return
                end
                ewt[i + 1] = 1.0 / ewt[i + 1]
            end
        end
        tolsf = ETA * vmnorm(n, yh[1], ewt)
        if tolsf > 0.01
            tolsf = tolsf * 200.0
            if nst == 0
                fprintf(__stderrp, "lsoda -- at start of problem, too much accuracy\n")
                fprintf(__stderrp, "         requested for precision of machine,\n")
                fprintf(__stderrp, "         suggested scaling factor = %g\n", tolsf)
                terminate(istate)
                freevectors()
                return
            end
            fprintf(__stderrp, "lsoda -- at t = %g, too much accuracy requested\n", t[])
            fprintf(__stderrp, "         for precision of machine, suggested\n")
            fprintf(__stderrp, "         scaling factor = %g\n", tolsf)
            istate[] = -2
            terminate2(y, t)
            return
        end
        if tn + h == tn
            nhnil += 1
            if nhnil <= mxhnil
                fprintf(__stderrp, "lsoda -- warning..internal t = %g and h = %g are\n", tn, h)
                fprintf(__stderrp, "         such that in the machine, t + h = t on the next step\n")
                fprintf(__stderrp, "         solver will continue anyway.\n")
                if nhnil == mxhnil
                    fprintf(__stderrp, "lsoda -- above warning has been issued %d times,\n", nhnil)
                    fprintf(__stderrp, "         it will not be issued again for this problem\n")
                end
            end
        end
        stoda(neq, y, f, _data)
        if kflag == 0
            init = 1
            if meth != mused
                tsw = tn
                maxord = mxordn
                if meth == 2
                    maxord = mxords
                end
                jstart = -1
                if ixpr
                    if meth == 2
                        fprintf(__stderrp, "[lsoda] a switch to the stiff method has occurred ")
                    end
                    if meth == 1
                        fprintf(__stderrp, "[lsoda] a switch to the nonstiff method has occurred")
                    end
                    fprintf(__stderrp, "at t = %g, tentative step size h = %g, step nst = %d\n", tn, h, nst)
                end
            end
            if itask == 1
                if (tn - tout) * h < 0.0
                    continue
                end
                @c intdy(tout, 0, y, &iflag)
                t[] = tout
                istate[] = 2
                illin = 0
                freevectors()
                return
            end
            if itask == 2
                successreturn(y, t, itask, ihit, tcrit, istate)
                return
            end
            if itask == 3
                if (tn - tout) * h >= 0.0
                    successreturn(y, t, itask, ihit, tcrit, istate)
                    return
                end
                continue
            end
            if itask == 4
                if (tn - tout) * h >= 0.0
                    @c intdy(tout, 0, y, &iflag)
                    t[] = tout
                    istate[] = 2
                    illin = 0
                    freevectors()
                    return
                else
                    hmx = fabs(tn) + fabs(h)
                    ihit = fabs(tn - tcrit) <= (100.0ETA) * hmx
                    if ihit
                        successreturn(y, t, itask, ihit, tcrit, istate)
                        return
                    end
                    tnext = tn + h * (1.0 + 4.0ETA)
                    if (tnext - tcrit) * h <= 0.0
                        continue
                    end
                    h = (tcrit - tn) * (1.0 - 4.0ETA)
                    jstart = -2
                    continue
                end
            end
            if itask == 5
                hmx = fabs(tn) + fabs(h)
                ihit = fabs(tn - tcrit) <= (100.0ETA) * hmx
                successreturn(y, t, itask, ihit, tcrit, istate)
                return
            end
        end
        if kflag == -1 || kflag == -2
            fprintf(__stderrp, "lsoda -- at t = %g and step size h = %g, the\n", tn, h)
            if kflag == -1
                fprintf(__stderrp, "         error test failed repeatedly or\n")
                fprintf(__stderrp, "         with fabs(h) = hmin\n")
                istate[] = -4
            end
            if kflag == -2
                fprintf(__stderrp, "         corrector convergence failed repeatedly or\n")
                fprintf(__stderrp, "         with fabs(h) = hmin\n")
                istate[] = -5
            end
            big = 0.0
            imxer = 1
            begin
                size = fabs(acor[i + 1]) * ewt[i + 1]
                if big < size
                    big = size
                    imxer = i
                end
            end
            terminate2(y, t)
            return
        end
    end
end
function stoda(neq, y, f, _data)
    local corflag
    local orderflag
    local i
    local i1
    local j
    local m
    local ncf
    local del
    local delp
    local dsm
    local dup
    local exup
    local r
    local rh
    local rhup
    local told
    local pdh
    local pnorm
    kflag = 0
    told = tn
    ncf = 0
    ierpj = 0
    iersl = 0
    jcur = 0
    delp = 0.0
    if jstart == 0
        lmax = maxord + 1
        nq = 1
        l = 2
        ialth = 2
        rmax = 10000.0
        rc = 0.0
        el0 = 1.0
        crate = 0.7
        hold = h
        nslp = 0
        ipup = miter
        icount = 20
        irflag = 0
        pdest = 0.0
        pdlast = 0.0
        ratio = 5.0
        cfode(2)
        cm2[i + 1] = (tesco[i + 1])[2] * (elco[i + 1])[(i + 1) + 1]
        cfode(1)
        cm1[i + 1] = (tesco[i + 1])[2] * (elco[i + 1])[(i + 1) + 1]
        resetcoeff()
    end
    if jstart == -1
        ipup = miter
        lmax = maxord + 1
        if ialth == 1
            ialth = 2
        end
        if meth != mused
            cfode(meth)
            ialth = l
            resetcoeff()
        end
        if h != hold
            rh = h / hold
            h = hold
            @c scaleh(&rh, &pdh)
        end
    end
    if jstart == -2
        if h != hold
            rh = h / hold
            h = hold
            @c scaleh(&rh, &pdh)
        end
    end
    while 1
        while 1
            if fabs(rc - 1.0) > ccmax
                ipup = miter
            end
            if nst >= nslp + msbp
                ipup = miter
            end
            tn += h
            begin
                yp1 = yh[i1 + 1]
                yp2 = yh[(i1 + 1) + 1]
                yp1[i + 1] += yp2[i + 1]
            end
            pnorm = vmnorm(n, yh[1], ewt)
            @c correction(neq, y, f, &corflag, pnorm, &del, &delp, &told, &ncf, &rh, &m, _data)
            if corflag == 0
                break
            end
            if corflag == 1
                rh = @warn("max(rh,hmin/fabs(h))")
                @c scaleh(&rh, &pdh)
                continue
            end
            if corflag == 2
                kflag = -2
                hold = h
                jstart = 1
                return
            end
        end
        jcur = 0
        if m == 0
            dsm = del / (tesco[nq + 1])[2]
        end
        if m > 0
            dsm = vmnorm(n, acor, ewt) / (tesco[nq + 1])[2]
        end
        if dsm <= 1.0
            kflag = 0
            nst += 1
            hu = h
            nqu = nq
            mused = meth
            begin
                yp1 = yh[j + 1]
                r = el[j + 1]
                yp1[i + 1] += r * acor[i + 1]
            end
            icount -= 1
            if icount < 0
                @c methodswitch(dsm, pnorm, &pdh, &rh)
                if meth != mused
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                    rmax = 10.0
                    endstoda()
                    break
                end
            end
            ialth -= 1
            if ialth == 0
                rhup = 0.0
                if l != lmax
                    yp1 = yh[lmax + 1]
                    savf[i + 1] = acor[i + 1] - yp1[i + 1]
                    dup = vmnorm(n, savf, ewt) / (tesco[nq + 1])[3]
                    exup = 1.0 / Cdouble(l + 1)
                    rhup = 1.0 / (1.4 * pow(dup, exup) + 1.4e-6)
                end
                @c orderswitch(&rhup, dsm, &pdh, &rh, &orderflag)
                if orderflag == 0
                    endstoda()
                    break
                end
                if orderflag == 1
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                    rmax = 10.0
                    endstoda()
                    break
                end
                if orderflag == 2
                    resetcoeff()
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                    rmax = 10.0
                    endstoda()
                    break
                end
            end
            if ialth > 1 || l == lmax
                endstoda()
                break
            end
            yp1 = yh[lmax + 1]
            yp1[i + 1] = acor[i + 1]
            endstoda()
            break
        else
            kflag -= 1
            tn = told
            begin
                yp1 = yh[i1 + 1]
                yp2 = yh[(i1 + 1) + 1]
                yp1[i + 1] -= yp2[i + 1]
            end
            rmax = 2.0
            if fabs(h) <= hmin * 1.00001
                kflag = -1
                hold = h
                jstart = 1
                break
            end
            if kflag > -3
                rhup = 0.0
                @c orderswitch(&rhup, dsm, &pdh, &rh, &orderflag)
                if orderflag == 1 || orderflag == 0
                    if orderflag == 0
                        rh = @warn("min(rh,0.2)")
                    end
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                end
                if orderflag == 2
                    resetcoeff()
                    rh = @warn("max(rh,hmin/fabs(h))")
                    @c scaleh(&rh, &pdh)
                end
                continue
            else
                if kflag == -10
                    kflag = -1
                    hold = h
                    jstart = 1
                    break
                else
                    rh = 0.1
                    rh = @warn("max(hmin/fabs(h),rh)")
                    h *= rh
                    yp1 = yh[1]
                    y[i + 1] = yp1[i + 1]
                    (tn, y + 1, savf + 1, _data)
                    nfe += 1
                    yp1 = yh[2]
                    yp1[i + 1] = h * savf[i + 1]
                    ipup = miter
                    ialth = 5
                    if nq == 1
                        continue
                    end
                    nq = 1
                    l = 2
                    resetcoeff()
                    continue
                end
            end
        end
    end
end
function ewset(itol, rtol, atol, ycur)
    local i
    @switch itol begin
            @case 1
            ewt[i + 1] = rtol[1] * fabs(ycur[i + 1]) + atol[1]
            break
            @case 2
            ewt[i + 1] = rtol[1] * fabs(ycur[i + 1]) + atol[i + 1]
            break
            @case 3
            ewt[i + 1] = rtol[i + 1] * fabs(ycur[i + 1]) + atol[1]
            break
            @case 4
            ewt[i + 1] = rtol[i + 1] * fabs(ycur[i + 1]) + atol[i + 1]
            break
        end
end
function intdy(t, k, dky, iflag)
    local i
    local ic
    local j
    local jj
    local jp1
    local c
    local r
    local s
    local tp
    iflag[] = 0
    if k < 0 || k > nq
        fprintf(__stderrp, "[intdy] k = %d illegal\n", k)
        iflag[] = -1
        return
    end
    tp = (tn - hu) - (100.0ETA) * (tn + hu)
    if (t - tp) * (t - tn) > 0.0
        fprintf(__stderrp, "intdy -- t = %g illegal. t not in interval tcur - hu to tcur\n", t)
        iflag[] = -2
        return
    end
    s = (t - tn) / h
    ic = 1
    ic *= jj
    c = Cdouble(ic)
    yp1 = yh[l + 1]
    dky[i + 1] = c * yp1[i + 1]
    begin
        jp1 = j + 1
        ic = 1
        ic *= jj
        c = Cdouble(ic)
        yp1 = yh[jp1 + 1]
        dky[i + 1] = c * yp1[i + 1] + s * dky[i + 1]
    end
    if k == 0
        return
    end
    r = pow(h, Cdouble(-k))
    dky[i + 1] *= r
end
function cfode(meth)
    local i
    local nq
    local nqm1
    local nqp1
    local agamq
    local fnq
    local fnqm1
    pc = 13
    local pint
    local ragq
    local rqfac
    local rq1fac
    local tsign
    local xpin
    if meth == 1
        (elco[1])[1] = 1.0
        (elco[1])[2] = 1.0
        (tesco[1])[1] = 0.0
        (tesco[1])[2] = 2.0
        (tesco[2])[1] = 1.0
        (tesco[12])[3] = 0.0
        pc[1] = 1.0
        rqfac = 1.0
        begin
            rq1fac = rqfac
            rqfac = rqfac / Cdouble(nq)
            nqm1 = nq - 1
            fnqm1 = Cdouble(nqm1)
            nqp1 = nq + 1
            pc[nq + 1] = 0.0
            pc[i + 1] = pc[(i - 1) + 1] + fnqm1 * pc[i + 1]
            pc[1] = fnqm1 * pc[1]
            pint = pc[1]
            xpin = pc[1] / 2.0
            tsign = 1.0
            begin
                tsign = -tsign
                pint += (tsign * pc[i + 1]) / Cdouble(i)
                xpin += (tsign * pc[i + 1]) / Cdouble(i + 1)
            end
            (elco[nq + 1])[1] = pint * rq1fac
            (elco[nq + 1])[2] = 1.0
            (elco[nq + 1])[(i + 1) + 1] = (rq1fac * pc[i + 1]) / Cdouble(i)
            agamq = rqfac * xpin
            ragq = 1.0 / agamq
            (tesco[nq + 1])[2] = ragq
            if nq < 12
                (tesco[nqp1 + 1])[1] = (ragq * rqfac) / Cdouble(nqp1)
            end
            (tesco[nqm1 + 1])[3] = ragq
        end
        return
    end
    pc[1] = 1.0
    rq1fac = 1.0
    begin
        fnq = Cdouble(nq)
        nqp1 = nq + 1
        pc[nqp1 + 1] = 0.0
        pc[i + 1] = pc[(i - 1) + 1] + fnq * pc[i + 1]
        pc[1] *= fnq
        (elco[nq + 1])[i + 1] = pc[i + 1] / pc[2]
        (elco[nq + 1])[2] = 1.0
        (tesco[nq + 1])[1] = rq1fac
        (tesco[nq + 1])[2] = Cdouble(nqp1) / (elco[nq + 1])[1]
        (tesco[nq + 1])[3] = Cdouble(nq + 2) / (elco[nq + 1])[1]
        rq1fac /= fnq
    end
    return
end
function scaleh(rh, pdh)
    local r
    local j
    local i
    rh[] = @warn("min(*rh,rmax)")
    rh[] = rh[] / @warn("max(1.,fabs(h)*hmxi**rh)")
    if meth == 1
        irflag = 0
        pdh[] = @warn("max(fabs(h)*pdlast,0.000001)")
        if (rh[] * pdh[]) * 1.00001 >= sm1[nq + 1]
            rh[] = sm1[nq + 1] / pdh[]
            irflag = 1
        end
    end
    r = 1.0
    begin
        r *= rh[]
        yp1 = yh[j + 1]
        yp1[i + 1] *= r
    end
    h *= rh[]
    rc *= rh[]
    ialth = l
end
function prja(neq, y, f, _data)
    local i
    local ier
    local j
    local fac
    local hl0
    local r
    local r0
    local yj
    nje += 1
    ierpj = 0
    jcur = 1
    hl0 = h * el0
    if miter != 2
        fprintf(__stderrp, "[prja] miter != 2\n")
        return
    end
    if miter == 2
        fac = vmnorm(n, savf, ewt)
        r0 = (((1000.0 * fabs(h)) * ETA) * Cdouble(n)) * fac
        if r0 == 0.0
            r0 = 1.0
        end
        begin
            yj = y[j + 1]
            r = @warn("max(sqrteta*fabs(yj),r0/ewt[j])")
            y[j + 1] += r
            fac = -hl0 / r
            (tn, y + 1, acor + 1, _data)
            (wm[i + 1])[j + 1] = (acor[i + 1] - savf[i + 1]) * fac
            y[j + 1] = yj
        end
        nfe += n
        pdnorm = fnorm(n, wm, ewt) / fabs(hl0)
        (wm[i + 1])[i + 1] += 1.0
        @c dgefa(wm, n, ipvt, &ier)
        if ier != 0
            ierpj = 1
        end
        return
    end
end
function vmnorm(n, v, w)
    local i
    local vm
    vm = 0.0
    vm = @warn("max(vm,fabs(v[i])*w[i])")
    return vm
end
function fnorm(n, a, w)
    local i
    local j
    local an
    local sum
    local ap1
    an = 0.0
    begin
        sum = 0.0
        ap1 = a[i + 1]
        sum += fabs(ap1[j + 1]) / w[j + 1]
        an = @warn("max(an,sum*w[i])")
    end
    return an
end
function correction(neq, y, f, corflag, pnorm, del, delp, told, ncf, rh, m, _data)
    local i
    local rm
    local rate
    local dcon
    m[] = 0
    corflag[] = 0
    rate = 0.0
    del[] = 0.0
    yp1 = yh[1]
    y[i + 1] = yp1[i + 1]
    (tn, y + 1, savf + 1, _data)
    nfe += 1
    while 1
        if m[] == 0
            if ipup > 0
                prja(neq, y, f, _data)
                ipup = 0
                rc = 1.0
                nslp = nst
                crate = 0.7
                if ierpj != 0
                    corfailure(told, rh, ncf, corflag)
                    return
                end
            end
            acor[i + 1] = 0.0
        end
        if miter == 0
            yp1 = yh[2]
            begin
                savf[i + 1] = h * savf[i + 1] - yp1[i + 1]
                y[i + 1] = savf[i + 1] - acor[i + 1]
            end
            del[] = vmnorm(n, y, ewt)
            yp1 = yh[1]
            begin
                y[i + 1] = yp1[i + 1] + el[1] * savf[i + 1]
                acor[i + 1] = savf[i + 1]
            end
        else
            yp1 = yh[2]
            y[i + 1] = h * savf[i + 1] - (yp1[i + 1] + acor[i + 1])
            solsy(y)
            del[] = vmnorm(n, y, ewt)
            yp1 = yh[1]
            begin
                acor[i + 1] += y[i + 1]
                y[i + 1] = yp1[i + 1] + el[1] * acor[i + 1]
            end
        end
        if del[] <= (100.0pnorm) * ETA
            break
        end
        if m[] != 0 || meth != 1
            if m[] != 0
                rm = 1024.0
                if del[] <= 1024.0 * delp[]
                    rm = del[] / delp[]
                end
                rate = @warn("max(rate,rm)")
                crate = @warn("max(0.2*crate,rm)")
            end
            dcon = (del[] * @warn("min(1.,1.5*crate)")) / ((tesco[nq + 1])[2] * conit)
            if dcon <= 1.0
                pdest = @warn("max(pdest,rate/fabs(h*el[1]))")
                if pdest != 0.0
                    pdlast = pdest
                end
                break
            end
        end
        (m[])[]
        if m[] == maxcor || m[] >= 2 && del[] > 2.0 * delp[]
            if miter == 0 || jcur == 1
                corfailure(told, rh, ncf, corflag)
                return
            end
            ipup = miter
            m[] = 0
            rate = 0.0
            del[] = 0.0
            yp1 = yh[1]
            y[i + 1] = yp1[i + 1]
            (tn, y + 1, savf + 1, _data)
            nfe += 1
        else
            delp[] = del[]
            (tn, y + 1, savf + 1, _data)
            nfe += 1
        end
    end
end
function corfailure(told, rh, ncf, corflag)
    local j
    local i1
    local i
    ncf += 1
    rmax = 2.0
    tn = told[]
    begin
        yp1 = yh[i1 + 1]
        yp2 = yh[(i1 + 1) + 1]
        yp1[i + 1] -= yp2[i + 1]
    end
    if fabs(h) <= hmin * 1.00001 || ncf[] == mxncf
        corflag[] = 2
        return
    end
    corflag[] = 1
    rh[] = 0.25
    ipup = miter
end
function solsy(y)
    iersl = 0
    if miter != 2
        printf("solsy -- miter != 2\n")
        return
    end
    if miter == 2
        dgesl(wm, n, ipvt, y, 0)
    end
    return
end
function methodswitch(dsm, pnorm, pdh, rh)
    local lm1
    local lm1p1
    local lm2
    local lm2p1
    local nqm1
    local nqm2
    local rh1
    local rh2
    local rh1it
    local exm2
    local dm2
    local exm1
    local dm1
    local alpha
    local exsm
    if meth == 1
        if nq > 5
            return
        end
        if dsm <= (100.0pnorm) * ETA || pdest == 0.0
            if irflag == 0
                return
            end
            rh2 = 2.0
            nqm2 = @warn("min(nq,mxords)")
        else
            exsm = 1.0 / Cdouble(l)
            rh1 = 1.0 / (1.2 * pow(dsm, exsm) + 1.2e-6)
            rh1it = 2.0rh1
            pdh[] = pdlast * fabs(h)
            if pdh[] * rh1 > 1.0e-5
                rh1it = sm1[nq + 1] / pdh[]
            end
            rh1 = @warn("min(rh1,rh1it)")
            if nq > mxords
                nqm2 = mxords
                lm2 = mxords + 1
                exm2 = 1.0 / Cdouble(lm2)
                lm2p1 = lm2 + 1
                dm2 = vmnorm(n, yh[lm2p1 + 1], ewt) / cm2[mxords + 1]
                rh2 = 1.0 / (1.2 * pow(dm2, exm2) + 1.2e-6)
            else
                dm2 = dsm * (cm1[nq + 1] / cm2[nq + 1])
                rh2 = 1.0 / (1.2 * pow(dm2, exsm) + 1.2e-6)
                nqm2 = nq
            end
            if rh2 < ratio * rh1
                return
            end
        end
        rh[] = rh2
        icount = 20
        meth = 2
        miter = jtyp
        pdlast = 0.0
        nq = nqm2
        l = nq + 1
        return
    end
    exsm = 1.0 / Cdouble(l)
    if mxordn < nq
        nqm1 = mxordn
        lm1 = mxordn + 1
        exm1 = 1.0 / Cdouble(lm1)
        lm1p1 = lm1 + 1
        dm1 = vmnorm(n, yh[lm1p1 + 1], ewt) / cm1[mxordn + 1]
        rh1 = 1.0 / (1.2 * pow(dm1, exm1) + 1.2e-6)
    else
        dm1 = dsm * (cm2[nq + 1] / cm1[nq + 1])
        rh1 = 1.0 / (1.2 * pow(dm1, exsm) + 1.2e-6)
        nqm1 = nq
        exm1 = exsm
    end
    rh1it = 2.0rh1
    pdh[] = pdnorm * fabs(h)
    if pdh[] * rh1 > 1.0e-5
        rh1it = sm1[nqm1 + 1] / pdh[]
    end
    rh1 = @warn("min(rh1,rh1it)")
    rh2 = 1.0 / (1.2 * pow(dsm, exsm) + 1.2e-6)
    if rh1 * ratio < 5.0rh2
        return
    end
    alpha = @warn("max(0.001,rh1)")
    dm1 *= pow(alpha, exm1)
    if dm1 <= (1000.0ETA) * pnorm
        return
    end
    rh[] = rh1
    icount = 20
    meth = 1
    miter = 0
    pdlast = 0.0
    nq = nqm1
    l = nq + 1
end
function endstoda()
    local r
    local i
    r = 1.0 / (tesco[nqu + 1])[2]
    acor[i + 1] *= r
    hold = h
    jstart = 1
end
function orderswitch(rhup, dsm, pdh, rh, orderflag)
    local newq
    local i
    local exsm
    local rhdn
    local rhsm
    local ddn
    local exdn
    local r
    orderflag[] = 0
    exsm = 1.0 / Cdouble(l)
    rhsm = 1.0 / (1.2 * pow(dsm, exsm) + 1.2e-6)
    rhdn = 0.0
    if nq != 1
        ddn = vmnorm(n, yh[l + 1], ewt) / (tesco[nq + 1])[1]
        exdn = 1.0 / Cdouble(nq)
        rhdn = 1.0 / (1.3 * pow(ddn, exdn) + 1.3e-6)
    end
    if meth == 1
        pdh[] = @warn("max(fabs(h)*pdlast,0.000001)")
        if l < lmax
            rhup[] = @warn("min(*rhup,sm1[l]/*pdh)")
        end
        rhsm = @warn("min(rhsm,sm1[nq]/*pdh)")
        if nq > 1
            rhdn = @warn("min(rhdn,sm1[nq-1]/*pdh)")
        end
        pdest = 0.0
    end
    if rhsm >= rhup[]
        if rhsm >= rhdn
            newq = nq
            rh[] = rhsm
        else
            newq = nq - 1
            rh[] = rhdn
            if kflag < 0 && rh[] > 1.0
                rh[] = 1.0
            end
        end
    else
        if rhup[] <= rhdn
            newq = nq - 1
            rh[] = rhdn
            if kflag < 0 && rh[] > 1.0
                rh[] = 1.0
            end
        else
            rh[] = rhup[]
            if rh[] >= 1.1
                r = el[l + 1] / Cdouble(l)
                nq = l
                l = nq + 1
                yp1 = yh[l + 1]
                yp1[i + 1] = acor[i + 1] * r
                orderflag[] = 2
                return
            else
                ialth = 3
                return
            end
        end
    end
    if meth == 1
        if (rh[] * pdh[]) * 1.00001 < sm1[newq + 1]
            elseif kflag == 0 && rh[] < 1.1
                ialth = 3
                return
            end
        end
    else
        if kflag == 0 && rh[] < 1.1
            ialth = 3
            return
        end
    end
    if kflag <= -2
        rh[] = @warn("min(*rh,0.2)")
    end
    if newq == nq
        orderflag[] = 1
        return
    end
    nq = newq
    l = nq + 1
    orderflag[] = 2
end
function resetcoeff()
    local i
    local ep1
    ep1 = elco[nq + 1]
    el[i + 1] = ep1[i + 1]
    rc = (rc * el[1]) / el0
    el0 = el[1]
    conit = 0.5 / Cdouble(nq + 2)
end
function freevectors()
end
function _freevectors()
    local i
    if wm
        free(wm[i + 1])
    end
    if yh
        free(yh[i + 1])
    end
    free(yh)
    free(wm)
    free(ewt)
    free(savf)
    free(acor)
    free(ipvt)
    g_nyh = (g_lenyh = 0)
    yh = 0
    wm = 0
    ewt = 0
    savf = 0
    acor = 0
    ipvt = 0
end
function n_lsoda(y, n, x, xout, eps, yscal, devis, data)
    local i
    local istate
    local itask
    local _y
    local atol
    local rtol
    _y = Ptr{Cdouble}(calloc(3 * (n + 1), nothing))
    atol = (_y + n) + 1
    rtol = (atol + n) + 1
    begin
        _y[i + 1] = y[(i - 1) + 1]
        atol[i + 1] = eps * yscal[(i - 1) + 1]
    end
    istate = if init
            2
        else
            1
        end
    itask = 2
    @c lsoda(devis, n, _y, x, xout, 2, rtol, atol, itask, &istate, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, data)
    y[(i - 1) + 1] = _y[i + 1]
    free(_y)
    return if istate < 0
            istate
        else
            0
        end
end
function n_lsoda_terminate()
    if init
        _freevectors()
    end
    init = 0
end
