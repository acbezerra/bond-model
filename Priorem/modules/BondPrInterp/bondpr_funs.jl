function get_cvm_bond_price(svm, Vt, ttm, sigma; vb=nothing)
	# ####################################
	# ######## Extract Parameters ########
	# ####################################
	# Capital Structure
	m = svm.m
	c = svm.c
	p = svm.p
	
	# Default
	alpha = svm.pm.alpha
	pii = svm.pm.pi
	
	# Dynamics
	r = svm.pm.r
	gross_delta = svm.pm.gross_delta
	
	# Liquidity
	xi = svm.pm.xi
	kappa = svm.pm.kappa
	# ####################################
	
	if vb==nothing
		vb = get_cvm_vb(svm, sigma)
	end

	return cvm_bond_price(vb, log(Vt/vb), ttm, m, c, p, sigma, r,
	                        gross_delta, xi, kappa, alpha, pii)

end


# ##########################################################
# #################### Fin Diff Method #####################
# ##########################################################
# 1. Notice the value of ttm passed to
# the function call must match the 
# value of ttm used to calibrate the 
# functions f0, f1, f2 and f3.
# 2. Notice also the order of the 
# arguments in the f2 and f3 functions:
# first vt, then vbhl.
# 3. Finally, vbhl in the f2 and f3 
# functions below is the ratio of the 
# post- to -pre-volatility shock
# bankruptcy barriers.

# ##########################################################
# ############# Prices at Different Maturities #############
# ##########################################################


# ##########################################################
# ##########################################################

function get_interp_values(vals, f11_sitp, f12_sitp, f13_sitp,
                           f21_sitp, f22_sitp; 
		           fd_method=true)

	if fd_method
	    return f11_sitp[vals["vt"]],
		   f12_sitp[vals["vt"], vals["vbhl"]], 
	           f13_sitp[vals["vt"], vals["vbhl"]],
                   f21_sitp[vals["vt"]], f22_sitp[vals["vt"]]
	else
	    return f11_sitp[vals["ttm"]], f12_sitp[vals["ttm"]],
                   f13_sitp[vals["ttm"]], f21_sitp[vals["ttm"]],
                   f22_sitp[vals["ttm"]]
	end
end

### WARNING: ttm must equal m if fd_method=true!!!!
function bondpr(svm, Vt, vbl, ttm, vmax,
        	f11_sitp, f12_sitp, f13_sitp,
                f21_sitp, f22_sitp; c=nothing, p=nothing,
		fd_method=true)
	
        # ####################################
        # ######## Extract Parameters ########
        # ####################################
        # Capital Structure
        m = svm.m

    
        if c==nothing
            c = svm.c
        end
        if p==nothing
            p = svm.p
        end
    
        # Default
        alpha = svm.pm.alpha
        pii = svm.pm.pi
    
        # Dynamics
        r = svm.pm.r
        gross_delta = svm.pm.gross_delta
        iota = svm.pm.iota
    
        # Liquidity
        xi = svm.pm.xi
        kappa = svm.pm.kappa
    
        # Volatility
        _lambda = svm.pm.lambda
        sigmal = svm.pm.sigmal
        sigmah = svm.pm.sigmah
        # ####################################
    
	vt = log(Vt/float(vbl))
   
        # Default Barrier
        vbh = zhi_vb(m, c, p, sigmah, r,
                        gross_delta, iota, xi, kappa, alpha, pii)
    
        if vt < 0
                return on_default_payoff(vt, vbl, ttm,
                                         m, c, p, r, xi,
                                         kappa, alpha)
        else
            if vt > vmax
                   return rfbond_price(ttm, c, p, r, xi, kappa)
            else
    
    		if fd_method
    			vals = Dict("vt" => vt, 
    				    "vbhl" => vbh/float(vbl))
    		else
    			vals = Dict("ttm" => ttm)
    		end
    
    		f11, f12, f13, f21, f22 = get_interp_values(vals, f11_sitp,
    			         		            f12_sitp, f13_sitp,
                                                            f21_sitp, f22_sitp,
    			         		            fd_method=fd_method)
    
                # Maturity or Default prior to Volatility Shock:
                cf0 = no_vol_shock_cf_pv(vt, vbl, ttm,
                                         m, c, p, sigmal,
                                         r, gross_delta,
                                         xi, kappa, alpha, _lambda)
 
                # Volatility Shock Prior to Maturity:
                cf1 = c/rdisc(r, xi, kappa) * f11 +
                      (p - c/rdisc(r, xi, kappa)) * f12 +
                      (alpha * vbh/float(m) - c/rdisc(r, xi, kappa)) * f13

                cf2 = c/rdisc(r, xi, kappa) * f21 +
                      (p - c/rdisc(r, xi, kappa)) * f22 
 
                return min(cf0 + cf1 + cf2, 
		           rfbond_price(ttm, c, p, r, xi, kappa))
            end
        end
end



# Objective is to find V such that the value of the
# newly-issued (tau = m) risky bond price when sigma = sigmah
# is sufficiently close to the credit-risk-free bond price.
function get_bond_Vmax(svm; initV=nothing, tol=.5 * 1e-3, print=false)

	# ####################################
	# ######## Extract Parameters ########
	# ####################################
	# Capital Structure
	m = svm.m
	c = svm.c
	p = svm.p
	
	# Default
	alpha = svm.pm.alpha
	pii = svm.pm.pi
	
	# Dynamics
	r = svm.pm.r
	gross_delta = svm.pm.gross_delta
	iota = svm.pm.iota
	
	# Liquidity
	xi = svm.pm.xi
	kappa = svm.pm.kappa
	
	# Volatility
	sigmal = svm.pm.sigmal
	sigmah = svm.pm.sigmah
	# ####################################


	if initV == nothing
		initV = svm.pm.V0
	end

	bondVmax = 1.25 * initV 
	vb = zhi_vb(m, c, p, sigmah, r, gross_delta, iota, xi, kappa, alpha, pii)
	
	vmax = log(bondVmax/vb)
	rfbond = rfbond_price(m, c, p, r, xi, kappa)
	bondpr = cvm_bond_price(vb, vmax, m, m, c, p, sigmah, r,
	                            gross_delta, xi, kappa, alpha, pii)
	
	per_diff = (rfbond - bondpr) / rfbond
	
	cond = per_diff > tol 
	while cond
	    bondVmax = 1.025 * bondVmax
	    vmax = log(bondVmax / vb)
	    bondpr = cvm_bond_price(vb, vmax, m, m, c, p, sigmah, r,
	                        gross_delta, xi, kappa, alpha, pii)
	    per_diff = (rfbond - bondpr) / rfbond
	    cond = per_diff > tol
	end
	
	if print
	    println(string("Bond Percentage Difference: ", per_diff))
	end
	
	return bondVmax
end

