
function get_cvm_eq(V, sigma, svm) 
   v = log(V/get_cvm_vb(svm, sigma))
   return zhi_eq(v, svm.pm.m, svm.c, svm.p, 
                 sigma, svm.pm.r, svm.pm.gross_delta, 
                 svm.pm.iota, svm.pm.xi, svm.pm.kappa, 
                 svm.pm.alpha, svm.pm.pi)
end

function get_eq_Vmax(svm; initV=nothing, print_diffs=false)

	# ####################################
	# ######## Extract Parameters ########
	# ####################################
	# Capital Structure
	m = svm.pm.m
	c = svm.c
	C = get_param(svm, "C")
	p = svm.p
	
	# Default
	alpha = svm.pm.alpha
	pii = svm.pm.pi
	
	# Dynamics
	r = svm.pm.r
	delta = get_param(svm, "delta")
	# gross_delta = svm.pm.gross_delta
	# iota = svm.pm.iota
	
	# Liquidity
	xi = svm.pm.xi
	kappa = svm.pm.kappa
	
	# Volatility
	sigmah = svm.pm.sigmah
	# sigmal = svm.pm.sigmal
	# ####################################

	
	if initV == nothing
		initV = svm.V0
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
        println("Vmax: ", string(Vmax))

        # Limiting Equity Function:
	rfbond = rfbond_price(m, c, p, r, xi, kappa)

	phi0 = (- (1 - pii) * C + rfbond - p)/r
        println("phi0: ", string(phi0))

        phi1 = (delta/(r - get_rgrow(svm)))
        println("phi1: ", string(phi1))

        pv_rfdebt = get_pv_rfdebt(svm)
        println("pv_rfdebt: ", string(pv_rfdebt))

	debt_abs_per_diff = abs((get_cvm_debt_price(svm, sigmah; Vt=Vmax) -
                                    pv_rfdebt)/pv_rfdebt)
        eq_abs_per_diff = abs((get_cvm_eq(Vmax, sigmah, svm) -
                                 (phi0 + phi1 * Vmax))/(phi0 + phi1 * Vmax))

        cond = true
        count = 0
        while cond
            Vmax = 1.05*Vmax
            debt_abs_per_diff = abs((get_cvm_debt_price(svm, sigmah; Vt=Vmax) -
                                        pv_rfdebt) / pv_rfdebt)
            eq_abs_per_diff = abs((get_cvm_eq(Vmax, sigmah, svm) -
                                      (phi0 + phi1 * Vmax)) / (phi0 + phi1 * Vmax))
            cond = max(debt_abs_per_diff, eq_abs_per_diff) > 1e-3

            if print_diffs
                println("count: ", string(count), ", Vmax: ", string(Vmax))
                println("debt_abs_per_diff: ",  string(debt_abs_per_diff))
                println("eq_abs_per_diff: ", str(eq_abs_per_diff))
                println("cond: ",  string(cond))
	    end
	end

            count += 1

        println("debt_abs_per_diff: ", string(debt_abs_per_diff))
        println("eq_abs_per_diff:", string(eq_abs_per_diff))
        return Vmax

end
