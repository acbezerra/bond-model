
# ##################################################################
# ########################### f0_int ###############################
# ################################################################## 

function f0_mat(vt, vmax, _lambda, sigma, 
                c, p, r, gross_delta, xi, kappa;
		ttm_max=1.2, ttmN=20, uN=1e3)

	_, ttm_grid = grid_creator(0.0, ttm_max, ttmN)
	
	f0_int_vec = fetch(@spawn [f0_int(vt, vmax, _lambda, sigma, 
                                    c, p, r, gross_delta, xi, kappa; 
				    ttm=tau, N=uN) for tau=ttm_grid])
	f0_int_itp = interpolate(f0_int_vec, BSpline(Cubic(Line())), OnGrid())
	return Interpolations.scale(f0_int_itp, ttm_grid)
end


# ##################################################################
# ########################### f1_int ###############################
# ################################################################## 

function f1_mat(vt, _lambda, sigmal,
		r, gross_delta, xi, kappa;
		vmax=1.2, vN=1e3, ttm_max=1.2, uN=1e3)
		

    dv, vgrid = grid_creator(0.0, vmax, vN)
    du, ugrid = grid_creator(0.0, ttm_max, uN) 

    _, f1v_int_vec = f1v_int(vt,  _lambda, sigmal,
                       	     r, gross_delta, xi, kappa;
                                 vgrid=vgrid, ugrid=ugrid)
 
    f1v_itp = interpolate(f1v_int_vec, BSpline(Cubic(Line())), OnGrid())
    f1v_sitp = Interpolations.scale(f1v_itp, ugrid)


    f1v_int_u = fetch(@spawn [integral_u(f1v_sitp, u, zN=10^4) for u=ugrid])
    f1v_itp2 = interpolate(f1v_int_u, BSpline(Cubic(Line())), OnGrid())
    return Interpolations.scale(f1v_itp2, ugrid)
end


# ##################################################################
# ########################### f2_int ###############################
# ################################################################## 

function f2_mat(vt, vbhl, _lambda, sigmal, sigmah,
		r, gross_delta, xi, kappa;
		vmax=1.2, vN=1e3, ttm_max=1.2, ttmN=20, uN=1e3)
	
	dv, vgrid = grid_creator(0.0, vmax, vN)
	ttm_grid = linspace(0.0, ttm_max, ttmN)

	f2_surf = @spawn [f2_int(vt, vbhl, grid_creator(0.0, u, uN)[2],
				 vgrid, _lambda, sigmal, sigmah, r,
				 gross_delta, xi, kappa) for u=ttm_grid[2:end]]

	f2_itp = interpolate(fetch(f2_surf), BSpline(Cubic(Line())),
			    OnGrid())
	return Interpolations.scale(f2_itp, ttm_grid[2:end])
end


# ##################################################################
# ########################### f3_int ###############################
# ################################################################## 

function f3_mat(vt, vbhl, _lambda, sigmal, sigmah,
		r, gross_delta, xi, kappa;
		vmax=1.2, vN=1e3, ttm_max=1.2, ttmN=20, uN=1e3)
	
	
	dv, vgrid = grid_creator(0.0, vmax, vN)
	ttm_grid = linspace(0.0, ttm_max, ttmN)

	f3_surf = @spawn [f3_int(vt, vbhl, grid_creator(0.0, u, uN)[2], vgrid,
				 _lambda, sigmal, sigmah, r, gross_delta, xi,
				 kappa) for u=ttm_grid[2:end]]


	f3_itp = interpolate(fetch(f3_surf), BSpline(Cubic(Line())),
			    OnGrid())
	return Interpolations.scale(f3_itp, ttm_grid[2:end])
end


# ##################################################################
# ########################### Interp ###############################
# ################################################################## 

function mat_interp(vt, vbhl, svm;
       	vmax=1.2, vN=1e3, ttm_max=1.2, ttmN=20, uN=1e3)
	
        # ####################################
        # ######## Extract Parameters ########
        # ####################################
        # Capital Structure
        m = svm.pm.m
        c = svm.c
        p = svm.p
    
        # Taxes & Default
        
        # Dynamics
        r = svm.pm.r
        gross_delta = svm.pm.gross_delta
        # iota?
    
        # Liquidity
        xi = svm.pm.xi
        kappa = svm.pm.kappa
        
        # Volatility
        _lambda = svm.pm.lambda
        sigmal = svm.pm.sigmal
        sigmah = svm.pm.sigmah 
        # ####################################
        
	f0 = @spawn f0_mat(vt, vmax, _lambda, sigmal, 
			   c, p, r, gross_delta, xi,
			   kappa; ttm_max=ttm_max, ttmN=ttmN, uN=uN)

	f1 = @spawn f1_mat(vt, _lambda, sigmal, r, gross_delta, xi, kappa;
			   vmax=vmax, vN=vN, ttm_max=ttm_max, uN=uN)
		    
	f2 = @spawn f2_mat(vt, vbhl, _lambda, sigmal, sigmah, r, gross_delta,
			   xi, kappa; vmax=vmax, vN=vN, ttm_max=ttm_max,
			   ttmN=ttmN, uN=uN)
	
        f3 = @spawn f3_mat(vt, vbhl, _lambda, sigmal, sigmah, r, gross_delta,
			   xi, kappa; vmax=vmax, vN=vN, ttm_max=ttm_max,
			   ttmN=ttmN, uN=uN)
	
        return fetch(f0), fetch(f1), fetch(f2), fetch(f3)
end
