
function get_cvm_eq(svm, V, sigma;
                    mu_b::Float64=NaN, m::Float64=NaN,
                    c::Float64=NaN, p::Float64=NaN)
    mu_b, m, c, p = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=p)
    
    vb = get_cvm_vb(svm, sigma; mu_b=mu_b, m=m, c=c, p=p)

    return cvm_eq(log(V/float(vb)), mu_b, m, c, p, 
                  sigma, svm.pm.r, svm.pm.gross_delta, 
                  svm.pm.iota, svm.pm.xi, svm.pm.kappa, 
                  svm.pm.alpha, svm.pm.pi)
end


function get_eq_Vmax(svm; mu_b::Float64=NaN, m::Float64=NaN,
                     c::Float64=NaN, p::Float64=NaN,
                     initV::Float64=NaN, print_diffs::Bool=false)

    mu_b, m, c, p = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=p)
    
    if isnan(initV)      # use as vb
        initV = svm.pm.V0
    end

    # As shown in the Appendix, the slope of the
    # pre-volatility shock and constant volatility
    # equity functions converges to
    # delta/(r-rgrow)

    # As a first guess, set Vmax to the value of V
    # for which the probability of V crossing VB is
    # sufficiently close to zero:
    # Vmax = self.get_big_v_eq_max(vb=vb)
    Vmax = 1.25 * initV 
    println(string("Vmax: ", Vmax))

    # Limiting Equity Function:
    rfbond = rfbond_price(m, c, p, svm.pm.r,
                          svm.pm.xi, svm.pm.kappa)

    phi0 = mu_b * (- (1 - svm.pm.pi) * (m * c) + rfbond - p)/svm.pm.r
    println(string("phi0: ", phi0))

    phi1 = (get_param(svm, :delta)/(svm.pm.r - get_rgrow(svm)))
    println(string("phi1: ", phi1))

    pv_rfdebt = get_pv_rfdebt(svm, mu_b=mu_b, m=m, c=c, p=p)
    println("pv_rfdebt: ", string(pv_rfdebt))

    debt_abs_per_diff = abs((get_cvm_debt_price(svm, initV, svm.pm.sigmah;
                                                Vt=Vmax, mu_b=mu_b, m=m, c=c, p=p) -
                                pv_rfdebt)/pv_rfdebt)
    eq_abs_per_diff = abs((get_cvm_eq(svm, Vmax, svm.pm.sigmah;
                                      mu_b=mu_b, m=m, c=c, p=p) -
                             (phi0 + phi1 * Vmax))/(phi0 + phi1 * Vmax))

    cond = true
    count = 0
    while cond
        Vmax = 1.05 * Vmax
        debt_abs_per_diff = abs((get_cvm_debt_price(svm, initV, svm.pm.sigmah;
                                                    Vt=Vmax, mu_b=mu_b, m=m, c=c, p=p) -
                                    pv_rfdebt) / pv_rfdebt)
        eq_abs_per_diff = abs((get_cvm_eq(svm, Vmax, svm.pm.sigmah;
                                          mu_b=mu_b, m=m, c=c, p=p) -
                                  (phi0 + phi1 * Vmax)) / (phi0 + phi1 * Vmax))
        cond = maximum([debt_abs_per_diff, eq_abs_per_diff]) > 1e-3

        if print_diffs
            println(string("count: ", count, ", Vmax: ", Vmax))
            println(string("debt_abs_per_diff: ", debt_abs_per_diff))
            println(string("eq_abs_per_diff: ", eq_abs_per_diff))
            println(string("cond: ", cond))
	end
    end
    # ################
    
    count += 1

    println(string("debt_abs_per_diff: ", debt_abs_per_diff))
    println(string("eq_abs_per_diff: ", eq_abs_per_diff))
    
    return Vmax
end
