
# ##################################################################
# ######################### J2: f21, f22 ###########################
# ################################################################## 

function f21_mat(vt, vmax, _lambda, sigma, 
                r, gross_delta, xi, kappa;
		ttm_max=1.2, ttmN=20, uN=10^3)

	_, ttm_grid = grid_creator(0.0, ttm_max, ttmN)
	
	f21_int_vec = fetch(@spawn [f21_int(vt, vmax, _lambda, sigma, 
                                    r, gross_delta, xi, kappa; 
				    ttm=tau, N=uN) for tau=ttm_grid])
	f21_int_itp = interpolate(f21_int_vec, BSpline(Cubic(Line(OnGrid()))))
	return Interpolations.scale(f21_int_itp, ttm_grid)
end


function f22_mat(vt, vmax, _lambda, sigma, 
                r, gross_delta, xi, kappa;
		ttm_max=1.2, ttmN=20, uN=10^3)

	_, ttm_grid = grid_creator(0.0, ttm_max, ttmN)
	
	f22_int_vec = fetch(@spawn [f22_int(vt, vmax, _lambda, sigma, 
                                           r, gross_delta, xi, kappa; 
				           ttm=tau, N=uN) for tau=ttm_grid])
	f22_int_itp = interpolate(f22_int_vec, BSpline(Cubic(Line(OnGrid()))))
	return Interpolations.scale(f22_int_itp, ttm_grid)
end


# ##################################################################
# ########################### f1_int ###############################
# ################################################################## 

function f11_mat(vt, _lambda, sigmal,
		r, gross_delta, xi, kappa;
		vmax=1.2, vN=10^3, ttm_max=1.2, uN=10^3)
		

    dv, vgrid = grid_creator(0.0, vmax, vN)
    du, ugrid = grid_creator(0.0, ttm_max, uN) 

    _, f11v_int_vec = f11v_int(vt,  _lambda, sigmal,
                       	     r, gross_delta, xi, kappa;
                                 vgrid=vgrid, ugrid=ugrid)
 
    f11v_itp = interpolate(f11v_int_vec, BSpline(Cubic(Line(OnGrid()))))
    f11v_sitp = Interpolations.scale(f11v_itp, ugrid)


    f11v_int_u = fetch(@spawn [integral_u(f11v_sitp, u, zN=10^4) for u=ugrid])
    f11v_itp2 = interpolate(f11v_int_u, BSpline(Cubic(Line(OnGrid()))))
    return Interpolations.scale(f11v_itp2, ugrid)
end


# ##################################################################
# ########################### f2_int ###############################
# ################################################################## 

function f12_mat(vt, vbhl, _lambda, sigmal, sigmah,
		r, gross_delta, xi, kappa;
		vmax=1.2, vN=10^3, ttm_max=1.2, ttmN=20, uN=10^3)
	
	dv, vgrid = grid_creator(0.0, vmax, vN)
	ttm_grid = range(0.0, stop=ttm_max, length=ttmN)  # linspace(0.0, ttm_max, ttmN)

	f12_surf = @spawn [f12_int(vt, vbhl, grid_creator(0.0, u, uN)[2],
				 vgrid, _lambda, sigmal, sigmah, r,
				 gross_delta, xi, kappa) for u=ttm_grid[2:end]]

	f12_itp = interpolate(fetch(f12_surf), BSpline(Cubic(Line())),
			    OnGrid())
	return Interpolations.scale(f12_itp, ttm_grid[2:end])
end


# ##################################################################
# ########################### f3_int ###############################
# ################################################################## 

function f13_mat(vt, vbhl, _lambda, sigmal, sigmah,
		r, gross_delta, xi, kappa;
		vmax=1.2, vN=10^3, ttm_max=1.2, ttmN=20, uN=10^3)
	
	
	dv, vgrid = grid_creator(0.0, vmax, vN)
	ttm_grid = range(0.0, stop=ttm_max, length=ttmN) # linspace(0.0, ttm_max, ttmN)

	f13_surf = @spawn [f13_int(vt, vbhl, grid_creator(0.0, u, uN)[2], vgrid,
				 _lambda, sigmal, sigmah, r, gross_delta, xi,
				 kappa) for u=ttm_grid[2:end]]


	f13_itp = interpolate(fetch(f13_surf), BSpline(Cubic(Line())),
			    OnGrid())
	return Interpolations.scale(f13_itp, ttm_grid[2:end])
end


# ##################################################################
# ########################### Interp ###############################
# ################################################################## 

function mat_interp(vt, vbhl, svm;
       	vmax=1.2, vN=10^3, ttm_max=1.2, ttmN=20, uN=10^3)
	
        # ####################################
        # ######## Extract Parameters ########
        # ####################################
        # Capital Structure
        m = svm.m
    
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
        
	f11 = @spawn f11_mat(vt, _lambda, sigmal, r, gross_delta, xi, kappa;
			   vmax=vmax, vN=vN, ttm_max=ttm_max, uN=uN)
		    
	f12 = @spawn f12_mat(vt, vbhl, _lambda, sigmal, sigmah, r, gross_delta,
			   xi, kappa; vmax=vmax, vN=vN, ttm_max=ttm_max,
			   ttmN=ttmN, uN=uN)
	
        f13 = @spawn f13_mat(vt, vbhl, _lambda, sigmal, sigmah, r, gross_delta,
			   xi, kappa; vmax=vmax, vN=vN, ttm_max=ttm_max,
			   ttmN=ttmN, uN=uN)
    
	f21 = @spawn f21_mat(vt, vmax, _lambda, sigmal, 
			   r, gross_delta, xi,
			   kappa; ttm_max=ttm_max, ttmN=ttmN, uN=uN)
    
	f22 = @spawn f22_mat(vt, vmax, _lambda, sigmal, 
			   r, gross_delta, xi,
			   kappa; ttm_max=ttm_max, ttmN=ttmN, uN=uN)

	f11mat = fetch(f11)
	f12mat = fetch(f12)
	f13mat = fetch(f13)
	f21mat = fetch(f21)
	f22mat = fetch(f22)
    
	matfuns = Dict("f21"=> f21mat,
                       "f22"=> f22mat,
		       "f11"=> f11mat,
			"f12"=> f12mat,
			"f13"=> f13mat)
 
        return matfuns
end
