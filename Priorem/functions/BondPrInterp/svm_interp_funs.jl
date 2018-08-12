# SVM Interpolation Functions

# Functions from the Analytic Functions Module:
# rdisc
# rdisc_pvs
# rfbond_price
# psi_v_td
# dv_psi_v_td
# cvm_F
# cvm_G
# zhi_vb
# on_default_payoff
# no_vol_shock_cf_pv


function f0_int(vt, vmax, _lambda, sigma, 
                c, p, r, gross_delta, xi, kappa; 
		ttm=1.0, N=10e3, ugrid=nothing)

    if ugrid == nothing
        du =  (ttm - 1e-4)/N 
 	ugrid=linspace(du, ttm, N) 
    else
	du = ugrid[2] - ugrid[1]
    end

    f0_integrand = @spawn [_lambda * exp(-rdisc_pvs(r, xi, kappa, _lambda) * u) *
                            rfbond_price(ttm - u, c, p, r, xi, kappa) *
    			psi_v_td(vt, vmax, u, sigma, r, gross_delta) for u=ugrid]
    return sum(fetch(f0_integrand)) * du
end


function f1_int(vt, vgrid, _lambda, sigmal,
          	r, gross_delta, xi, kappa;
    	        ttm=1.0, N=10e3, ugrid=nothing)

    dv = vgrid[2] - vgrid[1] 
    if ugrid == nothing
        du =  (ttm - 1e-4)/N 
 	ugrid=linspace(du, ttm, N) 
    else
	du = ugrid[2] - ugrid[1]
    end

    f1_integrand = @spawn [_lambda * exp(- rdisc_pvs(r, xi, kappa, _lambda) * u) * 
			   (- dv_psi_v_td(vt, v, u, sigmal, r, gross_delta))
			   for u=ugrid, v=vgrid]
    return sum(fetch(f1_integrand)) * du * dv
end


function f2_int(vt, vbhl, ugrid, vgrid,
                _lambda, sigmal, sigmah,
    	    r, gross_delta, xi, kappa)

    # vbhl is the ratio vbh/vbl

    ttm = ugrid[end]
    du = ugrid[2] - ugrid[1] 
    dv = vgrid[2] - vgrid[1] 

    f2_integrand = @spawn [_lambda * exp(- rdisc_pvs(r, xi, kappa, _lambda) * u) *
                               (exp(-rdisc(r, xi, kappa) * (ttm - u)) * 
    			    (1 - cvm_F(v - log(vbhl), ttm - u, sigmah, r, gross_delta))) *
                               (-dv_psi_v_td(vt, v, u, sigmal, r, gross_delta)) *
                               (v - log(vbhl) > 0)
                              for u=ugrid, v=vgrid]
    return sum(fetch(f2_integrand)) * du * dv
end


function f3_int(vt, vbhl, ugrid, vgrid, #ttm, #  vmax,
                _lambda, sigmal, sigmah,
                r, gross_delta,
                xi, kappa)

    # vbhl is the ratio vbh/vbl

    ttm = ugrid[end]
    du = ugrid[2] - ugrid[1] 
    dv = vgrid[2] - vgrid[1] 

    f3_integrand = @spawn [_lambda * exp(- rdisc_pvs(r, xi, kappa, _lambda) * u) *
                               cvm_G(v - log(vbhl), ttm - u,
                                     sigmah, r, gross_delta, xi, kappa) *
                               (-dv_psi_v_td(vt, v, u, sigmal, r, gross_delta)) *
                               (v - log(vbhl) > 0)
                              for u=ugrid, v=vgrid]
    return sum(fetch(f3_integrand)) * du * dv
end


function bondpr_int_funs(vtgrid, vbhlgrid,
    		         ugrid, vgrid,
			 c, p,
                         _lambda, sigmal, sigmah,
                         r, gross_delta,
			 xi, kappa)

    # ttm = ugrid[end]
    dt = (ugrid[end] - ugrid[1])/length(ugrid)
    dv = (vgrid[end] - vgrid[1])/length(vgrid)


    vmax = vgrid[end]
    ttm = ugrid[end]
    N = 10^4
    dt2 = ttm/N
    ugrid2=linspace(dt2, ttm, N)
   
    # Function of vt alone, since ttm is given!
    f0_int_vec_future = @spawn [f0_int(vt, ugrid2, vmax,
       				   _lambda, sigmal, 
    		                   c, p, r, 
    				   gross_delta, xi, kappa)
    			   for vt=vtgrid]
    f0_int_vec = fetch(f0_int_vec_future)

   # Function of vt alone.
   f1_int_vec_future = @spawn [f1_int(vt, ugrid, vgrid, 
   	                           _lambda, sigmal, 
   				   r, gross_delta,
                                      xi, kappa)
   			    for vt=vtgrid]
   f1_int_vec = fetch(f1_int_vec_future)

   # Function of vt and vbhl:
   f2_int_surf_future = @spawn [f2_int(vt, vbhl,
   				    ugrid, vgrid,
                                       _lambda, sigmal, sigmah,
   	                            r, gross_delta, 
   				    xi, kappa)
   			    for vbhl=vbhlgrid, vt=vtgrid]
   f2_int_surf = fetch(f2_int_surf_future)


   # Function of vt and vbhl:
   f3_int_surf_future = @spawn [f3_int(vt, vbhl, ugrid, vgrid,
                                       _lambda, sigmal, sigmah,
                                       r, gross_delta,
                                       xi, kappa)
                                for vbhl=vbhlgrid, vt=vtgrid]
   f3_int_surf  = fetch(f3_int_surf_future)

    return f0_int_vec, f1_int_vec, f2_int_surf, f3_int_surf

end


function bondpr_funs_interpolator(vbhlgrid, vbhlgrid_refined_size,
    			      vtgrid, vtgrid_refined_size,
    			      f0_int_vec, f1_int_vec,
    			      f2_int_surf, f3_int_surf)
      

        vbhlgrid_refined = linspace(vbhlgrid[1], vbhlgrid[end], vbhlgrid_refined_size)
        vtgrid_refined = linspace(vtgrid[1], vtgrid[end], vtgrid_refined_size)
        
        # Interpolate the functions: 
        f0_itp = interpolate(f0_int_vec, BSpline(Cubic(Line())), OnGrid())
        f0_sitp = Interpolations.scale(f0_itp, vtgrid)
        
        f1_itp = interpolate(f1_int_vec, BSpline(Cubic(Line())), OnGrid())
        f1_sitp = Interpolations.scale(f1_itp, vtgrid) # Scale
        
        f2_itp = interpolate(f2_int_surf, BSpline(Cubic(Line())), OnGrid())
        f2_sitp = Interpolations.scale(f2_itp, vbhlgrid, vtgrid) # Scale
        
        f3_itp = interpolate(f3_int_surf, BSpline(Cubic(Line())), OnGrid())
        f3_sitp = Interpolations.scale(f3_itp, vbhlgrid, vtgrid) 

        return f0_sitp, f1_sitp, f2_sitp, f3_sitp
end


function bondpr_interp(svm, Vt, vbl, ttm, vmax,
    		   f0_sitp, f1_sitp, f2_sitp, f3_sitp)

    # ####################################
    # ######## Extract Parameters ########
    # ####################################
    # Capital Structure
    m = svm.pm.m
    c = svm.c
    p = svm.p

    # Default
    alpha = svm.pm.alpha
    pi = svm.pm.pi

    # Dynamics
    r = svm.pm.r
    gross_delta = svm.pm.gross_delta
    iota = svm.pm.iota

    # Liquidity
    xi = svm.pm.xi
    k = svm.pm.kappa

    # Volatility
    _lambda = svm.pm.lambda
    sigmal = svm.pm.sigmal
    sigmah = svm.pm.sigmah
    # ####################################

    vt = log(Vt/vbl)

    # Default Barrier
    vbh = zhi_vb(m, c, p, sigmah, r,
                    gross_delta, iota, xi, k, alpha, pi)

    if vt < 0
            return on_default_payoff(vt, vbl, ttm,
                                     m, c, p, r, xi,
                                     k, alpha)
    else
        if vt > vmax
            return rfbond_price(ttm, c, p, r, xi, k)
        else
            # Maturity or Default prior to Volatility Shock:
            cf0 = no_vol_shock_cf_pv(vt, vbl, ttm,
                                     m, c, p, sigmal,
                                     r, gross_delta,
                                     xi, k, alpha, _lambda)

            # Volatility Shock Prior to Maturity:
            cf1 = c/rdisc(r, xi, k) * f1_sitp[vt] +
                  (p - c/rdisc(r, xi, k)) * f2_sitp[vbh/vbl, vt] +
                  (alpha * vbh/m - c/rdisc(r, xi, k)) * f3_sitp[vbh/vbl, vt]

          return min(cf0 + cf1 + f0_sitp[vt], rfbond_price(ttm, c, p, r, xi, k))
        end
    end
end

